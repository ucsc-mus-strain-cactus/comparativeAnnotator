####################################################################################################
# Configuration
# Modify variables below as new releases are made. Ideally, the only thing you should have to change
# is the MSCA_VERSION variable and the genomes if there are new genomes.
####################################################################################################

# Genomes in this analysis
#genomes = Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ CAROLIEiJ PAHARIEiJ
genomes = 129S1
refGenome = C57B6J

# Data directories
#MSCA_PROJ_DIR = /hive/groups/recon/projs/mus_strain_cactus
MSCA_PROJ_DIR = .
MSCA_DATA_DIR = ${MSCA_PROJ_DIR}/pipeline_data
MSCA_VERSION = 1412v2
GENCODE_VERSION = VM2
ANNOTATION_DIR = ${MSCA_DATA_DIR}/comparative/${MSCA_VERSION}/comparativeAnnotation
TRANS_MAP_DIR = ${MSCA_DATA_DIR}/comparative/${MSCA_VERSION}/transMap
ASSEMBLY_HUB_DIR = ${MSCA_DATA_DIR}/comparative/${MSCA_VERSION}/assemblyHub

# Input files
HAL = ${MSCA_DATA_DIR}/comparative/${MSCA_VERSION}/cactus/${MSCA_VERSION}.hal
# TODO: write a target for extracting this on the fly
ATTRS = /hive/users/ifiddes/mouse_release_data/wgEncodeGencodeAttrsVM2.attrs
#ATTRS = ${MSCA_PROJ_DIR}/data/wgEncodeGencodeBasicVM2.attrs

# jobTree configuration
batchSystem = parasol
maxThreads = 30
maxCpus = 1024
defaultMemory = 8589934592
maxJobDuration = 36000

# halTools configuration
h5prefix=~

# Set shell, environmental variables
SHELL = /bin/bash -e
export SHELLOPTS := pipefail
export PYTHONPATH := ./:${PYTHONPATH}
export PATH := ./bin/:./sonLib/bin:./submodules/jobTree/bin:./hal/bin/:${PATH}

# Create temporary file extension
host = $(shell hostname)
ppid = $(shell echo $$PPID)
tmpExt = ${host}.${ppid}.tmp

# Can run hgSql commands on Kolossus
ifneq ($(wildcard ${HOME}/.hg.rem.conf),)
    export HGDB_CONF=${HOME}/.hg.rem.conf
endif

# Options to pass to pslCDnaFilter
filterOpts = -localNearBest=0.0001

# Gencode src data files
srcGencodeSet = wgEncodeGencodeBasic${GENCODE_VERSION}
srcDataDir = ${TRANS_MAP_DIR}/data
srcBasicPre = ${srcDataDir}/${srcGencodeSet}
srcBasicGp = ${srcBasicPre}.gp
srcBasicCds = ${srcBasicPre}.cds
srcBasicBed = ${srcBasicPre}.bed
srcBasicPsl = ${srcBasicPre}.psl

# mapping files
mappedDataDir = ${TRANS_MAP_DIR}/mapped
mappedRegionIdPsls = ${genomes:%=${mappedDataDir}/%.region.idpsl}
mappedBlockPsls = ${genomes:%=${mappedDataDir}/%.block.psl}

# chained
chainedDataDir = ${TRANS_MAP_DIR}/results/chained
chainedPsls = ${genomes:%=${chainedDataDir}/%.chained.psl}

# filtered
filteredDataDir = ${TRANS_MAP_DIR}/results/filtered
filteredPsls = ${genomes:%=${filteredDataDir}/%.filtered.psl}
filteredPslStats = ${genomes:%=${filteredDataDir}/%.filtered.psl.basestats}

# sequence files
targetSequenceDir = ${TRANS_MAP_DIR}/data/genomes
targetFastaFiles = ${genomes:%=${targetSequenceDir}/%.fa}
targetTwoBitFiles = ${genomes:%=${targetSequenceDir}/%.2bit}
targetChromSizes = ${genomes:%=${targetSequenceDir}/%.chrom.sizes}
queryFasta = ${targetSequenceDir}/${refGenome}.fa
queryTwoBit = ${targetSequenceDir}/${refGenome}.2bit
queryChromSizes = ${targetSequenceDir}/${refGenome}.chrom.sizes

# gene-check
geneCheckDir = ${TRANS_MAP_DIR}/results/geneCheck
geneCheckGps = ${genomes:%=${geneCheckDir}/%.gp}
geneCheckEvals = ${genomes:%=${geneCheckDir}/%.gene-check}
geneCheckEvalsBed = ${genomes:%=${geneCheckDir}/%.gene-check.bed}

# annotation
targetBedFiles = ${geneCheckEvalsBed}
jobTreeDir = .${MSCA_VERSION}_jobTree
log = ${ANNOTATION_DIR}/log.txt

# assemblyHub
halJobTreeDir = .${MSCA_VERSION}_halJobTree