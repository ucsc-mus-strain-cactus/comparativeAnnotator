####################################################################################################
# Configuration
# Modify variables below as new releases are made. Ideally, the only thing you should have to change
# is the MSCA_VERSION variable and the genomes if there are new genomes.
####################################################################################################

# Genomes in this analysis
genomes = C57B6J Rattus 129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ CAROLIEiJ PAHARIEiJ
refGenome = C57B6J
refGenomeSQLName = mm10

# Data directories
#MSCA_PROJ_DIR = /hive/groups/recon/projs/mus_strain_cactus
MSCA_PROJ_DIR = .
MSCA_DATA_DIR = ${MSCA_PROJ_DIR}/pipeline_data
MSCA_VERSION = 1412v3
GENCODE_VERSION = VM4
GENOMES_DIR = ${MSCA_DATA_DIR}/assemblies/${MSCA_VERSION}
ANNOTATION_DIR = ${MSCA_DATA_DIR}/comparative/${MSCA_VERSION}/comparativeAnnotation
TRANS_MAP_DIR = ${MSCA_DATA_DIR}/comparative/${MSCA_VERSION}/transMap
ASSEMBLY_HUB_DIR = ${MSCA_DATA_DIR}/comparative/${MSCA_VERSION}/assemblyHub
METRICS_DIR = ${MSCA_DATA_DIR}/comparative/${MSCA_VERSION}/metrics

# Input files
HAL = ${MSCA_DATA_DIR}/comparative/${MSCA_VERSION}/cactus/${MSCA_VERSION}_dev.hal

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
ifneq ($(wildcard ${HOME}/.hg.conf),)
    export HGDB_CONF=${HOME}/.hg.conf
endif

# Options to pass to pslCDnaFilter
#filterOpts = -localNearBest=0.0001

# sequence files
targetFastaFiles = ${genomes:%=${GENOMES_DIR}/%.fa}
targetTwoBitFiles = ${genomes:%=${GENOMES_DIR}/%.2bit}
targetChromSizes = ${genomes:%=${GENOMES_DIR}/%.chrom.sizes}
queryFasta = ${GENOMES_DIR}/${refGenome}.fa
queryTwoBit = ${GENOMES_DIR}/${refGenome}.2bit
queryChromSizes = ${GENOMES_DIR}/${refGenome}.chrom.sizes

# Gencode src data files
srcGencodeSet = wgEncodeGencodeBasic${GENCODE_VERSION}
srcPseudoGeneSet = wgEncodeGencodePseudoGene${GENCODE_VERSION}
srcDataDir = ${TRANS_MAP_DIR}/data
srcAttrs = ${srcDataDir}/wgEncodeGencodeAttrs${GENCODE_VERSION}
srcBasicGp = ${srcDataDir}/${srcGencodeSet}.gp
srcPseudoGp = ${srcDataDir}/${srcPseudoGeneSet}.gp
srcBasicCds = ${srcDataDir}/${srcGencodeSet}.cds
srcPseudoCds = ${srcDataDir}/${srcPseudoGeneSet}.cds
srcBasicPsl = ${srcDataDir}/${srcGencodeSet}.psl
srcPseudoPsl = ${srcDataDir}/${srcPseudoGeneSet}.psl
srcCombinedGp = ${srcDataDir}/${GENCODE_VERSION}.BasicPseudoCombined.gp
srcCombinedPsl = ${srcDataDir}/${GENCODE_VERSION}.BasicPseudoCombined.psl
srcCombinedCds = ${srcDataDir}/${GENCODE_VERSION}.BasicPseudoCombined.cds
srcCombinedCheck = ${srcDataDir}/${GENCODE_VERSION}.BasicPseudoCombined.gene-check
srcCombinedCheckBed = ${srcDataDir}/${GENCODE_VERSION}.BasicPseudoCombined.gene-check.bed

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

# gene-check
geneCheckDir = ${TRANS_MAP_DIR}/results/geneCheck
geneCheckGps = ${genomes:%=${geneCheckDir}/%.gp}
geneCheckEvals = ${genomes:%=${geneCheckDir}/%.gene-check}
geneCheckEvalsBed = ${genomes:%=${geneCheckDir}/%.gene-check.bed}
# for Dent's legacy code comparisons
geneCheckDetails = ${genomes:%=${geneCheckDir}/%.gene-check-details}
geneCheckDetailsBed = ${genomes:%=${geneCheckDir}/%.gene-check-details.bed}

# annotation
targetBedFiles = ${geneCheckEvalsBed}
jobTreeDir = .${MSCA_VERSION}_jobTree
log = ${MSCA_VERSION}_log.txt

# assemblyHub
halJobTreeDir = .${MSCA_VERSION}_halJobTree
halLog = ${MSCA_VERSION}_assembly_hub_log.txt

