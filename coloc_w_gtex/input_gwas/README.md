# input_gwas
for generating the finemapped gwas (frpm dap-g) input to fastenloc
# Script
1. `gwas_inputgen.py`
- for generating the GWAS input (fine-mapped results) for fastenloc
- can be applied to any files with snpid, csid and pip columns
- Please use `python3 gwas_inputgen.py -h` for input arguments. An example is in `crf_dapg_inputgen.sh`, which was used to generate the fastenloc input `crf_dapg.gz` (log file `crf_dapg_inputgen.log`).