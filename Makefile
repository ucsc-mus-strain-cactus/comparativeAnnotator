# modify config.mk to modify the pipeline
#include config.mk
#include config_1412v3.mk
include config_1411.mk

all: init srcData mapping chaining filtered extractFasta geneCheck annotation assemblyHub plots

init:
	# TODO: why does this error out?
	#git submodule update --init
	cd sonLib && make
	cd jobTree && make
	cd hal && make

####################################################################################################
# Retrieve src data. Uses hgSql and related Kent tools.
####################################################################################################
srcData: ${srcAttrs} ${srcBasicGp} ${srcPseudoGp} ${srcBasicPsl} ${srcPseudoPsl} ${srcBasicCds} ${srcPseudoCds} ${srcCombinedPsl} ${srcCombinedCheckBed} ${srcCombinedCheck}

${srcAttrs}:
	@mkdir -p $(dir $@)
	hgsql -Ne 'select geneId,geneName,geneType,transcriptId,transcriptType from $(notdir $@)' ${refGenomeSQLName} > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

# awk expression to edit chrom names in UCSC format.  Assumse all alts are version 1.
# chr1_GL456211_random, chrUn_GL456239
editUcscChrom = $$chromCol=="chrM"{$$chromCol="MT"} {$$chromCol = gensub("_random$$","", "g", $$chromCol);$$chromCol = gensub("^chr.*_([0-9A-Za-z]+)$$","\\1.1", "g", $$chromCol);  gsub("^chr","",$$chromCol); print $$0}
${srcBasicGp}:
	@mkdir -p $(dir $@)
	hgsql -Ne 'select * from ${srcGencodeSet}' ${refGenomeSQLName} | cut -f 2- | tawk -v chromCol=2 '${editUcscChrom}' > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${srcPseudoGp}:
	@mkdir -p $(dir $@)
	hgsql -Ne 'select * from ${srcPseudoGeneSet}' ${refGenomeSQLName} | cut -f 2- | tawk -v chromCol=2 '${editUcscChrom}' > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${srcBasicCds}: ${srcBasicPsl}

${srcPseudoCds}: ${srcPseudoPsl}

${srcBasicPsl}: ${srcBasicGp}
	@mkdir -p $(dir $@)
	genePredToFakePsl mm10 ${srcGencodeSet} stdout ${srcBasicCds} | tawk -v chromCol=14 '${editUcscChrom}' > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${srcPseudoPsl}: ${srcPseudoGp}
	@mkdir -p $(dir $@)
	genePredToFakePsl mm10 ${srcPseudoGeneSet} stdout ${srcPseudoCds} | tawk -v chromCol=14 '${editUcscChrom}' > $@.${tmpExt}
	mv -f $@.${tmpExt} $@	

${srcCombinedPsl}: ${srcbasicPsl} ${srcPseudoPsl}
	@mkdir -p $(dir $@)
	cat ${srcBasicPsl} ${srcPseudoPsl} | sort -k 10 > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${srcCombinedGp}: ${srcPseudoGp} ${srcBasicGp}
	@mkdir -p $(dir $@)
	cat ${srcPseudoGp} ${srcBasicGp} > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${srcCombinedCds}: ${srcBasicCds} ${srcPseudoCds}
	@mkdir -p $(dir $@)
	cat ${srcPseudoCds} ${srcBasicCds} > $@.${tmpExt}
	mv -f $@.${tmpExt} $@	

${srcCombinedCheck}: ${srcCombinedGp} ${queryTwoBit}
	@mkdir -p $(dir $@)
	sort -k2,2 -k 4,4n $< | gene-check --allow-non-coding --genome-seqs=${queryTwoBit} stdin $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${srcCombinedCheckBed}: ${srcCombinedCheck}
	@mkdir -p $(dir $@)
	genePredCheckToBed ${srcCombinedGp} ${srcCombinedCheck} $@.${tmpExt}
	mv -f $@.${tmpExt} $@	

####################################################################################################
# Mapping. Also uses hgSql and related Kent tools.
####################################################################################################
mapping: ${mappedRegionIdPsls} ${mappedBlockPsls}
${mappedDataDir}/%.region.idpsl: ${srcCombinedCheckBed}
	@mkdir -p $(dir $@)
	halLiftover --tab --outPSLWithName ${HAL} ${refGenome} ${srcCombinedCheckBed} $* $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${mappedDataDir}/%.block.psl: ${mappedDataDir}/%.region.idpsl
	@mkdir -p $(dir $@)
	pslMap -mapFileWithInQName ${srcCombinedPsl} $< $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####################################################################################################
# Chaining alignments. Uses Kent tools, but does not need hgSql access.
####################################################################################################
chaining: ${chainedPsls} 

${chainedDataDir}/%.chained.psl: ${mappedDataDir}/%.block.psl
	@mkdir -p $(dir $@)
	simpleChain -outPsl $< $@.${tmpExt}
	mv -f $@.${tmpExt} $@


####################################################################################################
# Filtering chained alignments.
####################################################################################################
filtered: ${filteredPsls} ${filteredPslStats}

${filteredDataDir}/%.filtered.psl: ${chainedDataDir}/%.chained.psl
	@mkdir -p $(dir $@)
	(pslCDnaFilter ${filterOpts} $< stdout | pslQueryUniq > $@.${tmpExt}) 2> /dev/null
	mv -f $@.${tmpExt} $@

${filteredDataDir}/%.filtered.psl.basestats: ${filteredDataDir}/%.filtered.psl
	@mkdir -p $(dir $@)
	generateBaseStats $< $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####################################################################################################
# Generating sequence files.
####################################################################################################
extractFasta: ${targetFastaFiles} ${targetTwoBitFiles} ${targetChromSizes} ${queryFasta} ${queryTwoBit} ${queryChromSizes}

${GENOMES_DIR}/%.fa:
	@mkdir -p $(dir $@)
	n="$(shell basename $@ | cut -d "." -f 1)" ;\
	hal2fasta ${HAL} $$n > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${GENOMES_DIR}/%.2bit: ${GENOMES_DIR}/%.fa
	@mkdir -p $(dir $@)
	faToTwoBit $< $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${GENOMES_DIR}/%.chrom.sizes:
	@mkdir -p $(dir $@)
	n="$(shell basename $@ | cut -d "." -f 1)" ;\
	halStats --chromSizes $$n ${HAL} > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${queryFasta}:
	@mkdir -p $(dir $@)
	n="$(shell basename $@ | cut -d "." -f 1)" ;\
	hal2fasta ${HAL} $$n > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${queryTwoBit}: ${queryFasta}
	@mkdir -p $(dir $@)
	faToTwoBit ${queryFasta} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${queryChromSizes}: ${queryTwoBit}
	@mkdir -p $(dir $@)
	twoBitInfo ${queryTwoBit} stdout | sort -k2rn > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####################################################################################################
# Running gene-check
####################################################################################################
geneCheck: ${geneCheckGps} ${geneCheckEvals} ${geneCheckEvalsBed}

${geneCheckDir}/%.gp: ${filteredDataDir}/%.filtered.psl ${srcCombinedCds}
	@mkdir -p $(dir $@)
	mrnaToGene -keepInvalid -quiet -genePredExt -ignoreUniqSuffix -insertMergeSize=0 -cdsFile=${srcCombinedCds} $< $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${geneCheckDir}/%.gene-check: ${geneCheckDir}/%.gp ${GENOMES_DIR}/%.2bit
	@mkdir -p $(dir $@)
	sort -k2,2 -k 4,4n $< | gene-check --allow-non-coding --genome-seqs=${GENOMES_DIR}/$*.2bit stdin ${geneCheckDir}/$*.gene-check.${tmpExt}
	mv -f ${geneCheckDir}/$*.gene-check.${tmpExt} ${geneCheckDir}/$*.gene-check

${geneCheckDir}/%.gene-check.bed: ${geneCheckDir}/%.gene-check
	@mkdir -p $(dir $@)
	genePredCheckToBed ${geneCheckDir}/$*.gp ${geneCheckDir}/$*.gene-check $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####################################################################################################
# Annotation pipeline. Going to ssh to ku to use the cluster if batch system is parasol.
####################################################################################################
annotation: ${ANNOTATION_DIR}/DONE

${ANNOTATION_DIR}/DONE: ${geneCheckEvalsBed} ${srcCombinedCheckBed}
	@mkdir -p $(dir $@)
	if [ -d ${jobTreeDir} ]; then rm -rf ${jobTreeDir}; fi
	if [ "${batchSystem}" = "parasol" ]; then \
		cwd="$(shell pwd)" ;\
		ssh ku -t "cd $$cwd && export PYTHONPATH=./ && \
		export PATH=./bin/:./sonLib/bin:./submodules/jobTree/bin:${PATH} && \
		python src/annotationPipeline.py --refGenome ${refGenome} --genomes ${genomes} --sizes ${targetChromSizes} \
		--psls ${filteredPsls} --beds ${targetBedFiles} --fastas ${targetFastaFiles} --refTwoBit ${queryTwoBit} \
		--annotationBed ${srcCombinedCheckBed} --batchSystem ${batchSystem} --gencodeAttributeMap ${srcAttrs} \
		--defaultMemory ${defaultMemory} --jobTree ${jobTreeDir} --maxJobDuration ${maxJobDuration} \
		--maxThreads ${maxThreads} --stats --outDir ${ANNOTATION_DIR} &> ${log}" ;\
	else \
		python src/annotationPipeline.py --refGenome ${refGenome} --genomes ${genomes} --sizes ${targetChromSizes} \
		--psls ${filteredPsls} --beds ${targetBedFiles} --fastas ${targetFastaFiles} --refTwoBit ${queryTwoBit} \
		--annotationBed ${srcCombinedCheckBed} --batchSystem ${batchSystem} --gencodeAttributeMap ${srcAttrs} \
		--defaultMemory ${defaultMemory} --jobTree ${jobTreeDir} --maxJobDuration ${maxJobDuration} \
		--maxThreads ${maxThreads} --stats --outDir ${ANNOTATION_DIR} &> ${log} ;\
	fi
	touch $@

####################################################################################################
# Building assemblyHub. Can't be run on ku due to weird issues with halLodExtract. SSH to kolossus.
# TODO: this creates all of the fastas/2bits from the hal again, unnecessarily.
####################################################################################################
assemblyHub: ${ASSEMBLY_HUB_DIR}/DONE

${ASSEMBLY_HUB_DIR}/DONE: ${ANNOTATION_DIR}/DONE
	if [ -d ${halJobTreeDir} ]; then rm -rf ${halJobTreeDir}; fi
	if [ -d ${ASSEMBLY_HUB_DIR} ]; then rm -rf ${ASSEMBLY_HUB_DIR}; mkdir ${ASSEMBLY_HUB_DIR}; fi
	cwd="$(shell pwd)" ;\
	bigBedDirs="$(shell /bin/ls -1d ${ANNOTATION_DIR}/bedfiles/* | paste -s -d ",")" ;\
	ssh kolossus.sdsc.edu -t "cd $$cwd && export PYTHONPATH=./ && export \
	PATH=./bin/:./sonLib/bin:./submodules/jobTree/bin:${PATH} && \
	python hal/assemblyHub/hal2assemblyHub.py ${HAL} ${ASSEMBLY_HUB_DIR} \
	--finalBigBedDirs $$bigBedDirs --maxThreads=${maxThreads} --batchSystem=singleMachine \
	--defaultMemory=${defaultMemory} --jobTree ${halJobTreeDir} \
	--maxJobDuration ${maxJobDuration} --stats --shortLabel ${MSCA_VERSION} \
	--longLabel ${MSCA_VERSION} --hub ${MSCA_VERSION} &>> ${halLog}"
	touch $@

####################################################################################################
# Generating some plots.
####################################################################################################
plots: ${METRICS_DIR}/DONE

${METRICS_DIR}/DONE: ${ANNOTATION_DIR}/DONE
	@mkdir -p $(dir $@)
	python scripts/coverage_identity_ok_plots.py --outDir ${METRICS_DIR} --genomes ${genomes} \
	--comparativeAnnotationDir ${ANNOTATION_DIR} --header ${MSCA_VERSION} --attrs ${srcAttrs} \
	--annotationBed ${srcCombinedCheckBed}
	touch $@