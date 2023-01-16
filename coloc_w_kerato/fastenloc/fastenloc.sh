#!/bin/sh
../../software/fastenloc/src/fastenloc \
    -eqtl ../../data/keratoconus/fastenloc/pheno2/dapg/fastenloc.kerato.annotation.vcf.gz \
    -gwas ../../data/keratoconus/fastenloc/pheno1/crf_fastenloc_input.tsv.gz \
    -prefix kerato.crf \
    -t kerato \
    -thread 4 \
    -total_variants 9904132