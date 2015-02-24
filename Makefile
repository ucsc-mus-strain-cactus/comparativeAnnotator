batchSystem = parasol
maxThreads = 30
maxCpus = 1024
defaultMemory = 8589934592
1411jobTree = jobTree_1411
1412jobTree = jobTree_1412
1411halJobTree = halJobTree_1411
1412halJobTree = halJobTree_1412
1411log = log_1411.txt
1412log = log_1412.txt
maxJobDuration = 36000
h5prefix = ~

export PYTHONPATH:=./:${PYTHONPATH}
export PATH:=./sonLib/bin:./submodules/jobTree/bin:./hal/bin/:${PATH}

1412genomes = C57B6J Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ CAROLIEiJ PAHARIEiJ
1411genomes = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ
refGenome = C57B6J

1411dataDir = /hive/users/ifiddes/mouse_release_data/1411
1412dataDir = /hive/users/ifiddes/mouse_release_data/1412
annotationBed = /hive/users/ifiddes/mouse_release_data/wgEncodeGencodeBasicVM2.gene-check.bed
gencodeAttributeMap = /hive/users/ifiddes/mouse_release_data/wgEncodeGencodeAttrsVM2.attrs
1411hal = /cluster/home/jcarmstr/public_html/mouseBrowser_1411/1411.hal
1412hal = /hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1412/cactus/1412.hal
1411trackHub = trackHub_1411/
1412trackHub = trackHub_1412/


all :
	git submodule update --init
	cd sonLib && make
	cd jobTree && make
	cd hal && make

1411 : all
	if [ -d ${1411jobTree} ]; then rm -rf ${1411jobTree}; fi
	python src/main.py --refGenome ${refGenome} --genomes ${1411genomes} --annotationBed ${annotationBed} \
	--dataDir ${1411dataDir} --gencodeAttributeMap ${gencodeAttributeMap} --outDir 1411_output/ \
	--maxThreads=${maxThreads} --batchSystem=${batchSystem} --defaultMemory=${defaultMemory} \
	--jobTree ${1411jobTree} --logLevel DEBUG --maxCpus ${maxCpus} --maxJobDuration ${maxJobDuration} \
	--stats &> ${1411log}
	if [ -d ${1411halJobTree} ]; then rm -rf ${1411halJobTree}; fi
	if [ -d {1411trackHub} ]; then rm -rf ${1411trackHub}; fi
	bigBedDirs = /bin/ls -1d 1411_output/bedfiles/* | paste -s -d ","
	python hal/assemblyHub/hal2assemblyHub.py ${1411hal} ${1411trackHub} --finalBigBedDirs ${bigBedDirs} \
	--maxThreads=${maxThreads} --batchSystem=${batchSystem} --defaultMemory=${defaultMemory} \
	--jobTree ${1411halJobTree} --logLevel DEBUG --maxCpus ${maxCpus} --maxJobDuration ${maxJobDuration} \
	--stats --shortLabel 1411 --longLabel 1411 --hub 1411 &> ${log}

1412 : all
	if [ -d ${1412jobTree} ]; then rm -rf ${1412jobTree}; fi
	python src/main.py --refGenome ${refGenome} --genomes ${1412genomes} --annotationBed ${annotationBed} \
	--dataDir ${1412dataDir} --gencodeAttributeMap ${gencodeAttributeMap} --outDir 1412_output/ \
	--maxThreads=${maxThreads} --batchSystem=${batchSystem} --defaultMemory=${defaultMemory} \
	--jobTree ${1412jobTree} --logLevel DEBUG --maxCpus ${maxCpus} --maxJobDuration ${maxJobDuration} \
	--stats &> ${1412log}
	if [ -d ${1412halJobTree} ]; then rm -rf ${1412halJobTree}; fi
	if [ -d {1412trackHub} ]; then rm -rf ${1412trackHub}; fi
	bigBedDirs = /bin/ls -1d output/bedfiles/* | paste -s -d ","
	python hal/assemblyHub/hal2assemblyHub.py ${1412hal} ${1412trackHub} --finalBigBedDirs ${bigBedDirs} \
	--maxThreads=${maxThreads} --batchSystem=${batchSystem} --defaultMemory=${defaultMemory} \
	--jobTree ${1412halJobTree} --logLevel DEBUG --maxCpus ${maxCpus} --maxJobDuration ${maxJobDuration} \
	--stats --shortLabel 1412 --longLabel 1412 --hub 1412 &> ${log}	
