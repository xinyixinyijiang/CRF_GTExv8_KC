# coloc_w_gtex

This directory contains the scripts to perform the colocalization between finemapped [CRF GWAS](https://www.nature.com/articles/s42003-020-01497-w) and [finemapped GTEx v8 cis e/sQTL data](https://zenodo.org/record/3517189#.YxtTFi8w35h)

## Scripts

1. `fastenloc.sh`: run colocalization analysis between CRF and GTEx v8 cis-e/sQTL signals in one tissue

2. `array_crf_gtex.sh`: run fastenloc.sh in EDDIE (UoE Cluster Computers) for 49 GTEx tissues using `array_crf_gtex_parameters.txt`. `qsub fastenloc_array.sh` was used to submit the batch jobs.

3. `res_format.sh` & `res_format.py`: reformat the raw results from fastenloc. The reformatted results are in the folder `./res_format`.