# dapg_workspace

This directory contains the scripts to fine-map [CRF GWAS](https://www.nature.com/articles/s42003-020-01497-w) using [DAP-G](https://github.com/xqwen/dap)

## Scripts
1. `./dapg_run.py` runs DAP-G using Z-statistics and LD Information. This script was used in bash script `dapg_run_95.sh` (ld control default 0.25, cs 0.95, no -ens) using files in the folder `dapg_run` as parameters.

2. `./tools.py` includes functions to summarize the results from DAP-G.

## Used software
- [DAP-G](https://github.com/xqwen/dap)

## Environment
`./tools.py`: Python v3.10.6, pandas v1.5.1
