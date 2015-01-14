batchSystem = parasol
maxThreads = 10
defaultMemory = 8589934592
jobTree = .jobTree
log = log.txt

export PYTHONPATH:=./:${PYTHONPATH}
export PATH:=./sonLib/bin:./submodules/jobTree/bin:${PATH}

genomes = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ
#genomes = Rattus
refGenome = C57B6J

annotationBed = /hive/groups/recon/projs/mus_strain_cactus/data/gene_check/wgEncodeGencodeBasicVM2.coding.gene-check.bed
dataDir = /hive/users/ifiddes/mus_strain_cactus/datafiles
gencodeAttributeMap = /cluster/home/markd/compbio/gencode/mus_strain_cactus/cactusMapCheck/experiments/2014-07-17.simpleChain/data/wgEncodeGencodeAttrsVM2.attrs

all :
	cd sonLib && make
	cd jobTree && make
	python lib/twobit/check_if_installed.py; if [ $$? == 3 ]; then python lib/twobit/setup_twobit.py build; python lib/twobit/setup_twobit.py install; fi

run : all
	if [ -d ${jobTree} ]; then rm -rf ${jobTree}; fi
	python src/main.py --refGenome ${refGenome} --genomes ${genomes} --annotationBed ${annotationBed} \
	--dataDir ${dataDir} --gencodeAttributeMap ${gencodeAttributeMap} \
	--maxThreads=${maxThreads} --batchSystem=${batchSystem} --defaultMemory=${defaultMemory} \
	--jobTree ${jobTree} --overwriteDb --logInfo &> ${log}
