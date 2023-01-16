for i in 4 6 29 36 54 56 60 61 69 77 83 90 97 98 100 103 112 115
do
python ../../../software/polyfun/extract_snpvar.py --sumstats ../../../data/keratoconus/polytfun/input/loci${i}_kerato.tsv.gz --out ./output/loci${i}_kerato_withvar.txt.gz
done