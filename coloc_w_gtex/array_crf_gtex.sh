#!/bin/bash
#$ -t 1-2
#$ -cwd
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=10:00:00
#$ -N array_crf_gtex

DATAFILE=array_crf_gtex_parameters.txt
TSULIST=../data/gtex/tissues.txt  #Path to the tissuelist file
NTSU=49 #number of tissues in TSULIST

FASTENLOC=../software/fastenloc/src/fastenloc #executable fastenloc
QTLTYPE=$(awk "NR==$SGE_TASK_ID" $DATAFILE | cut -d' ' -f1) # QTL TYPE - for naming the output files
GWAS=./input_gwas/crf_dapg.gz # GWAS input for fastenloc - can be generated from fastenloc_gwasinput_gen.py
FMSOFTWARE=$(awk "NR==$SGE_TASK_ID" $DATAFILE | cut -d' ' -f2) #Name of the finemapping software, for naming the output files
OUTPUTDIR=${QTLTYPE}_${FMSOFTWARE}
CORES=4
NVARIANTS=9904132 # Total number of variants in GWAS

for n in $(seq 1 $NTSU)
do
    TSU=`sed -n "${n}p" ${TSULIST}`
    QTLANNO=../data/gtex/${QTLTYPE}/fastenloc.${QTLTYPE}.annotation.${TSU}.vcf.gz  # QTL ANNO files for fastenloc
    PREFIX=${QTLTYPE}.${FMSOFTWARE}.${TSU}
    echo TSU
    echo 'command' is: bash fastenloc.sh $FASTENLOC $QTLANNO $GWAS $TSU $PREFIX $CORES $NVARIANTS $OUTPUTDIR
    bash fastenloc.sh $FASTENLOC $QTLANNO $GWAS $TSU $PREFIX $CORES $NVARIANTS $OUTPUTDIR
    echo 'done'
done