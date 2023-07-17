# CLPP

`ployfun/polyfun.sh`: calculates snp level priors for fine-mapping.

`FINEMAP/FINEMAP.sh`: peforms keratoconus GWAS fine-mapping at CRF-keratoconus GWAS overlapping loci using snps' priors from ``polyfun``.

`clpp_FINEMAPres.py`: summarizes keraoconus GWAS fine-mapping results from ``FINEMAP`` and performs CLPP calculation.

## Used Software

- [FINEMAP](http://www.christianbenner.com/) v1.4
- [PolyFUN](https://github.com/omerwe/polyfun)

## Environment

``clpp_FINEMAPres.py`` : Python v3.10.6, pandas v1.5.1