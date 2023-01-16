name=95
python dapg_run.py \
    -dapg ../software/dap/dap_src/dap-g \
    -loci 1 116 \
    -loci_name ./dapg_run/locinames.txt \
    -dapg_z ./dapg_run/zpaths.txt \
    -dapg_ld ./dapg_run/ldpaths.txt \
    -dapg_log True \
    -dapg_msize 5 \
    -dapg_t 2 \
    -dapg_ldcont 0.25 \
    -output ./dapg_output_${name} \
    -getcs_script ../software/dap/utility/get_credible_set.pl \
    -cs_prob 0.95 > dapg_run_${name}.log

python -c "import os; import pandas as pd; from tools import summary_dapg_res; summary_dapg_res('./dapg_output_${name}', 0.95)"