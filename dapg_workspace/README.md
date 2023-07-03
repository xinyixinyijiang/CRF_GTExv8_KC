# dapg_workspace

This folder is for performing [CRF GWAS](https://www.nature.com/articles/s42003-020-01497-w) finemapping using [DAP-G](https://github.com/xqwen/dap)

## 1. dapg_run.py
Python script for running dap-g using Z-statistics and LD Information. This script was used in bash script `dapg_run_95.sh` (ld control default 0.25, cs 0.95, no -ens, log file `dapg_run_95.log`) using files in the folder `dapg_run` as parameters.

For arguments, please run:
```bash
python3 dapg_run.py -h
```

## 2. tools.py
This script includes the function for summarizing results from DAP-G.