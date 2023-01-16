#!/bin/sh
FASTENLOC=$1 #executable fastenloc
QTLANNO=$2
GWAS=$3
TSU=$4
PREFIX=$5
CORES=$6
NVARIANTS=$7
OUTPUTDIR=$8 #Path to the tissuelist file

mkdir -p ${OUTPUTDIR}

${FASTENLOC} \
  -eqtl ${QTLANNO} \
  -gwas ${GWAS} \
  -t ${TSU} \
  -prefix  ${PREFIX}\
  -thread ${CORES} \
  -total_variants ${NVARIANTS}

mv *${PREFIX}* ./${OUTPUTDIR}
echo ${PREFIX} finished

