# coloc_w_gtex
This folder is for performing the colocalisation between finemapped [CRF GWAS](https://www.nature.com/articles/s42003-020-01497-w) and [finemapped GTEx v8 cis e/sQTL data](https://zenodo.org/record/3517189#.YxtTFi8w35h)

# Scripts
1. `fastenloc.sh`: running colocalization analysis between CRF and GTEx v8 cis-e/sQTL signals in one tissue
1. `array_crf_gtex.sh`: running fastenloc.sh in EDDIE (UoE Cluster Computers) for 49 GTEx tissues using array_crf_gtex_parameters.txt as input. Log files are in `./log`. Usage: `qsub fastenloc_array.sh`.

1. `res_format.sh` & `res_format.py`: reformatting the raw results from fastenloc. Results are in `./res_format`. Usage: `bash res_format.sh`.