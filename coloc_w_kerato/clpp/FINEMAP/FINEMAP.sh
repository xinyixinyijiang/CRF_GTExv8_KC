FINEMAP=/exports/igmm/eddie/QTLgroup/PERSONAL/Xinyi/PhDproject/software/FINEMAP/finemap_v1.4_x86_64/finemap_v1.4_x86_64
masterfile=./master
${FINEMAP} --sss --in-files ${masterfile}  --flip-beta --prior-snps --force-n-samples
