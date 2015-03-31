# modify config.mk to modify the pipeline
include config.mk

all: init srcData mapping chaining filtered extractFasta geneCheck annotation assemblyHub

init:
	# TODO: why does this error out?
	#git submodule update --init
	cd sonLib && make
	cd jobTree && make
	cd hal && make

####################################################################################################
# Retrieve src data. Uses hgSql and related Kent tools.
####################################################################################################
srcData: ${srcBasicGp} ${srcBasicBed} ${srcBasicPsl} ${srcBasicCds} ${srcAttrs}

# awk expression to edit chrom names in UCSC format.  Assumse all alts are version 1.
# chr1_GL456211_random, chrUn_GL456239
editUcscChrom = $$chromCol=="chrM"{$$chromCol="MT"} {$$chromCol = gensub("_random$$","", "g", $$chromCol);$$chromCol = gensub("^chr.*_([0-9A-Za-z]+)$$","\\1.1", "g", $$chromCol);  gsub("^chr","",$$chromCol); print $$0}
${srcBasicGp}:
	@mkdir -p $(dir $@)
	hgsql -Ne 'select * from ${srcGencodeSet}' ${refGenomeSQLName} | cut -f 2- | tawk -v chromCol=2 '${editUcscChrom}' >$@.${tmpExt}
	mv -f $@.${tmpExt} $@

${srcBasicCds}: ${srcBasicPsl}

${srcBasicBed}: ${srcBasicGp}
	@mkdir -p $(dir $@)
	genePredToBed $< $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${srcBasicPsl}: ${srcBasicGp}
	@mkdir -p $(dir $@)
	genePredToFakePsl mm10 ${srcGencodeSet} stdout ${srcBasicCds} | tawk -v chromCol=14 '${editUcscChrom}' >$@.${tmpExt}
	mv -f $@.${tmpExt} $@

${srcAttrs}:
	@mkdir -p $(dir $@)
	hgsql -Ne 'select geneId,geneName,geneType,transcriptId,transcriptType from $(notdir $@)' ${refGenomeSQLName} > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####################################################################################################
# Mapping. Also uses hgSql and related Kent tools.
####################################################################################################
mapping: ${mappedRegionIdPsls} ${mappedBlockPsls}
${mappedDataDir}/%.region.idpsl: ${srcBasicBed}
	@mkdir -p $(dir $@)
	halLiftover --tab --outPSLWithName ${HAL} ${refGenome} ${srcBasicBed} $* $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${mappedDataDir}/%.block.psl: ${mappedDataDir}/%.region.idpsl
	pslMap -mapFileWithInQName ${srcBasicPsl} $< $@.${tmpExt}
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
	(pslCDnaFilter ${filterOpts} $< stdout | pslQueryUniq >$@.${tmpExt}) 2> /dev/null
	mv -f $@.${tmpExt} $@

${filteredDataDir}/%.filtered.psl.basestats: ${filteredDataDir}/%.filtered.psl
	@mkdir -p $(dir $@)
	generateBaseStats $< $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####################################################################################################
# Generating sequence files.
####################################################################################################
extractFasta: ${targetFastaFiles} ${targetTwoBitFiles} ${targetChromSizes} ${queryFasta} ${queryTwoBit} ${queryChromSizes}

${targetSequenceDir}/%.fa:
	@mkdir -p $(dir $@)
	n="$(shell basename $@ | cut -d "." -f 1)" ;\
	hal2fasta ${HAL} $$n > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${targetSequenceDir}/%.2bit: ${targetSequenceDir}/%.fa
	@mkdir -p $(dir $@)
	faToTwoBit $< $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${targetSequenceDir}/%.chrom.sizes:
	n="$(shell basename $@ | cut -d "." -f 1)" ;\
	halStats --chromSizes $$n ${HAL} > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${queryFasta}:
	n="$(shell basename $@ | cut -d "." -f 1)" ;\
	hal2fasta ${HAL} $$n > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${queryTwoBit}: ${queryFasta}
	faToTwoBit ${queryFasta} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

${queryChromSizes}: ${queryTwoBit}
	twoBitInfo ${queryTwoBit} stdout | sort -k2rn > $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####################################################################################################
# Running gene-check
####################################################################################################
geneCheck: ${geneCheckGps} ${geneCheckEvals} ${geneCheckEvalsBed}

${geneCheckDir}/%.gp: ${filteredDataDir}/%.filtered.psl ${srcBasicCds}
	@mkdir -p $(dir $@)
	mrnaToGene -keepInvalid -quiet -genePredExt -ignoreUniqSuffix -insertMergeSize=0 -cdsFile=${srcBasicCds} $< stdout | tawk '$$6<$$7' >$@.${tmpExt}
	mv -f $@.${tmpExt} $@

# pattern rules only execute once for mutiple targets
${geneCheckDir}/%.gene-check: ${geneCheckDir}/%.gp ${targetSequenceDir}/%.2bit
	@mkdir -p ${geneCheckDir}
	sort -k2,2 -k 4,4n $< | gene-check --allow-non-coding --genome-seqs=${targetSequenceDir}/$*.2bit stdin ${geneCheckDir}/$*.gene-check.${tmpExt}
	mv -f ${geneCheckDir}/$*.gene-check.${tmpExt} ${geneCheckDir}/$*.gene-check

${geneCheckDir}/%.gene-check.bed: ${geneCheckDir}/%.gene-check
	@mkdir -p $(dir $@)
	genePredCheckToBed ${geneCheckDir}/$*.gp ${geneCheckDir}/$*.gene-check $@.${tmpExt}
	mv -f $@.${tmpExt} $@

####################################################################################################
# Annotation pipeline. Going to ssh to ku to use the cluster if batch system is parasol.
####################################################################################################
annotation: ${ANNOTATION_DIR}/DONE

${ANNOTATION_DIR}/DONE: ${geneCheckEvalsBed}
	if [ -d ${jobTreeDir} ]; then rm -rf ${jobTreeDir}; fi
	if [ ! -d ${ANNOTATION_DIR} ]; then mkdir ${ANNOTATION_DIR}; fi
	if [ "${batchSystem}" = "parasol" ]; then \
		cwd="$(shell pwd)" ;\
		ssh ku -t "cd $$cwd && export PYTHONPATH=./ && export \
		PATH=./bin/:./sonLib/bin:./submodules/jobTree/bin:${PATH} && \
		python src/annotationPipeline.py --refGenome ${refGenome} --genomes ${genomes} \
		--psls ${filteredPsls} --beds ${targetBedFiles} --fastas ${targetFastaFiles} --refTwoBit ${queryTwoBit} \
		--annotationBed ${srcBasicBed} --batchSystem ${batchSystem} --gencodeAttributeMap ${srcAttrs} \
		--defaultMemory ${defaultMemory} --jobTree ${jobTreeDir} --maxJobDuration ${maxJobDuration} \
		--maxCpus ${maxCpus} --stats --outDir ${ANNOTATION_DIR} --sizes ${targetChromSizes} \
		--psls ${filteredPsls} --beds ${targetBedFiles} &> ${log}" ;\
	else \
		python src/annotationPipeline.py --refGenome ${refGenome} --genomes ${genomes} \
		--psls ${filteredPsls} --beds ${targetBedFiles} --fastas ${targetFastaFiles} --refTwoBit ${queryTwoBit} \
		--annotationBed ${srcBasicBed} --batchSystem ${batchSystem} --gencodeAttributeMap ${srcAttrs} \
		--defaultMemory ${defaultMemory} --jobTree ${jobTreeDir} --maxJobDuration ${maxJobDuration} \
		--maxThreads ${maxThreads} --stats --outDir ${ANNOTATION_DIR} --sizes ${targetChromSizes} \
		--psls ${filteredPsls} --beds ${targetBedFiles} &> ${log} ;\
	fi
	touch ${ANNOTATION_DIR}/DONE

####################################################################################################
# Building assemblyHub. Can't be run on ku due to weird issues with halLodExtract.
# TODO: this creates all of the fastas/2bits from the hal again, unnecessarily.
####################################################################################################
assemblyHub: ${ASSEMBLY_HUB_DIR}/DONE

${ASSEMBLY_HUB_DIR}/DONE: ${ANNOTATION_DIR}/DONE
	if [ -d ${halJobTreeDir} ]; then rm -rf ${halJobTreeDir}; fi
	bigBedDirs="$(shell /bin/ls -1d ${ANNOTATION_DIR}/bedfiles/* | paste -s -d ",")" ;\
	python hal/assemblyHub/hal2assemblyHub.py ${HAL} ${ASSEMBLY_HUB_DIR} \
	--finalBigBedDirs $$bigBedDirs --maxThreads=${maxThreads} --batchSystem=singleMachine \
	--defaultMemory=${defaultMemory} --jobTree ${halJobTreeDir} \
	--maxJobDuration ${maxJobDuration} --stats --shortLabel ${MSCA_VERSION} \
	--longLabel ${MSCA_VERSION} --hub ${MSCA_VERSION} &>> ${log}
	touch ${ASSEMBLY_HUB_DIR}/DONE