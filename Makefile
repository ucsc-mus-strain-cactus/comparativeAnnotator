batchSystem = singleMachine
maxThreads = 30
maxCpus = 1024
defaultMemory = 8589934592
jobTree = .jobTree
log = log.txt
maxJobDuration = 36000

export PYTHONPATH:=./:${PYTHONPATH}
export PATH:=./sonLib/bin:./submodules/jobTree/bin:hal/bin/:${PATH}
export h5prefix=-prefix=~

genomes = FVBNJ
#genomes = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ
refGenome = C57B6J

rootDir := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

dataDir = ${rootDir}/datafiles
annotationBed = ${dataDir}/wgEncodeGencodeBasicVM2.gene-check.bed
gencodeAttributeMap = ${dataDir}/wgEncodeGencodeAttrsVM2.attrs
hal = /cluster/home/jcarmstr/public_html/mouseBrowser_1411/1411.hal
trackHub = trackHub/
bedFiles = output/bedfiles

all :
	cd sonLib && make
	cd jobTree && make
	cd hal && make
	python lib/twobit/check_if_installed.py; if [ $$? == 3 ]; then python lib/twobit/setup_twobit.py build; python lib/twobit/setup_twobit.py install; fi

run : all
	if [ -d ${jobTree} ]; then rm -rf ${jobTree}; fi
	python src/main.py --refGenome ${refGenome} --genomes ${genomes} --annotationBed ${annotationBed} \
	--dataDir ${dataDir} --gencodeAttributeMap ${gencodeAttributeMap} \
	--maxThreads=${maxThreads} --batchSystem=${batchSystem} --defaultMemory=${defaultMemory} \
	--jobTree ${jobTree} --logLevel DEBUG --maxCpus ${maxCpus} --maxJobDuration ${maxJobDuration} \
	--stats &> ${log}
	python hal/assemblyHub/hal2assemblyHub.py ${hal} ${trackHub} --bedDirs ${bedFiles} \
	--maxThreads=${maxThreads} --batchSystem=${batchSystem} --defaultMemory=${defaultMemory} \
	--jobTree .halJobTree --logLevel DEBUG --maxCpus ${maxCpus} --maxJobDuration ${maxJobDuration} \
	--stats &> ${log}
