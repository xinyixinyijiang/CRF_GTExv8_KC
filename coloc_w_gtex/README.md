# coloc_w_gtex

This directory contains the scripts to perform the colocalization between fine-mapped [CRF GWAS](https://www.nature.com/articles/s42003-020-01497-w) and [fine-mapped GTEx v8 cis e/sQTL data](https://zenodo.org/record/3517189#.YxtTFi8w35h).

## Scripts

1. `fastenloc.sh` runs colocalization analysis between CRF and GTEx v8 cis-e/sQTL signals in one tissue. The input CRF GWAS summary statistics (`./input/gwas/crf_dapg.gz`) are in the folder `./input_gwas/`.

<!-- 2. `array_crf_gtex.sh`: run fastenloc.sh in EDDIE (UoE Cluster Computers) for 49 GTEx tissues using `array_crf_gtex_parameters.txt`. -->

2. `res_format.sh` & `res_format.py` reformat the raw results from ``fastenloc``. The reformatted results are in the folder `./res_format/`.

## Used Software 

- [fastENLOC](https://github.com/xqwen/fastenloc) v2.0

## Environment

``res_format.py`` : Python v3.10.6, pandas v1.5.1