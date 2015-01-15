batchSystem = parasol
maxThreads = 10
defaultMemory = 8589934592
jobTree = .jobTree
log = log.txt

export PYTHONPATH:=./:${PYTHONPATH}
export PATH:=./sonLib/bin:./submodules/jobTree/bin:${PATH}

genomes = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ
refGenome = C57B6J

rootDir := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

annotationBed = ${rootDir}/wgEncodeGencodeBasicVM2.coding.gene-check.bed
dataDir = ${rootDir}/datafiles
gencodeAttributeMap = ${rootDir}/wgEncodeGencodeAttrsVM2.attrs

all :
	cd sonLib && make
	cd jobTree && make
	python lib/twobit/check_if_installed.py; if [ $$? == 3 ]; then python lib/twobit/setup_twobit.py build; python lib/twobit/setup_twobit.py install; fi

run : all
	if [ -d ${jobTree} ]; then rm -rf ${jobTree}; fi
	python src/main.py --refGenome ${refGenome} --genomes ${genomes} --annotationBed ${annotationBed} \
	--dataDir ${dataDir} --gencodeAttributeMap ${gencodeAttributeMap} \
	--maxThreads=${maxThreads} --batchSystem=${batchSystem} --defaultMemory=${defaultMemory} \
	--jobTree ${jobTree} --overwriteDb &> ${log}

details : all
	python scripts/gene_check_details.py --refGenome ${refGenome} --genomes ${genomes} \
	--dataDir ${dataDir}
