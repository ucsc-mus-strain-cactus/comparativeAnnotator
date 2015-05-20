####################################################################################################
# Annotation pipeline. Going to ssh to ku to use the cluster if batch system is parasol.
# modify config.mk to modify the pipeline
####################################################################################################

include config_1411.mk

all: init annotation plots

init:
	# TODO: why does this error out?
	# git submodule update --init
	cd sonLib && make
	cd jobTree && make
	cd hal && make

annotation: ${ANNOTATION_DIR}/DONE

${ANNOTATION_DIR}/DONE: ${geneCheckGps} ${srcCombinedGp}
	@mkdir -p $(dir $@)
	if [ -d ${jobTreeDir} ]; then rm -rf ${jobTreeDir}; fi
	if [ "${batchSystem}" = "parasol" ]; then \
		cwd="$(shell pwd)" ;\
		ssh ku -t "cd $$cwd && export PYTHONPATH=./ && \
		export PATH=./bin/:./sonLib/bin:./jobTree/bin:${PATH} && \
		python src/annotationPipeline.py --refGenome ${refGenome} --genomes ${genomes} --sizes ${targetChromSizes} \
		--psls ${filteredPsls} --gps ${geneCheckGps} --fastas ${targetFastaFiles} --refTwoBit ${queryTwoBit} \
		--annotationGp ${srcCombinedGp} --batchSystem ${batchSystem} --gencodeAttributeMap ${srcAttrs} \
		--defaultMemory ${defaultMemory} --jobTree ${jobTreeDir} --maxJobDuration ${maxJobDuration} \
		--maxThreads ${maxThreads} --stats --outDir ${ANNOTATION_DIR} &> ${log}" ;\
	else \
		python src/annotationPipeline.py --refGenome ${refGenome} --genomes ${genomes} --sizes ${targetChromSizes} \
		--psls ${filteredPsls} --gps ${geneCheckGps} --fastas ${targetFastaFiles} --refTwoBit ${queryTwoBit} \
		--annotationGp ${srcCombinedGp} --batchSystem ${batchSystem} --gencodeAttributeMap ${srcAttrs} \
		--defaultMemory ${defaultMemory} --jobTree ${jobTreeDir} --maxJobDuration ${maxJobDuration} \
		--maxThreads ${maxThreads} --stats --outDir ${ANNOTATION_DIR} &> ${log} ;\
	fi
	touch $@

plots: ${METRICS_DIR}/DONE
# TODO: rewrite coverage_identity_ok_plots to handle genePreds
${METRICS_DIR}/DONE: ${ANNOTATION_DIR}/DONE
	@mkdir -p $(dir $@)
	python scripts/coverage_identity_ok_plots.py --outDir ${METRICS_DIR} --genomes ${genomes} \
	--comparativeAnnotationDir ${ANNOTATION_DIR} --header ${MSCA_VERSION} --attrs ${srcAttrs} \
	--annotationBed ${srcCombinedCheckBed}
	touch $@
