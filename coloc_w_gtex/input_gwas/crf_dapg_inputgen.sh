python gwas_inputgen.py \
    -gwas_fm_raw ../../dapg_workspace/dapg_output_95_pip_summary.tsv \
    -snpid_col snp \
    -locusid_col locus \
    -csid_col cs \
    -pip_col pip \
    -gtexid_ref ../../data/gtex/rsid_gtexsnpid_ref.gz \
    -output ./crf_dapg.gz > ./crf_dapg_inputgen.log