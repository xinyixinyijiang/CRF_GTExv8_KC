# CRF_GTEXV8_KC

This repo includes:

1. [coloc_w_gtex/](coloc_w_gtex/README.md) : colocalization between independent (i.e. fine-mapped) signals of corneal resistance factor (CRF) GWAS and cis-e/sQTL from GTEx v8.
2. [coloc_w_kerato/](coloc_w_kerato/README.md) : colocalization between independent signals of CRF GWAS and keratoconus GWAS. 
3. [dapg_workspace/](dapg_workspace/README.md) : fine-mapping of CRF GWAS and keratoconus GWAS using DAP-G
4. [scrnaseq_celltype_enrichment/](scrnaseq_celltype_enrichment/README.md) :  enrichment of CRF signals linked GTEx v8 e/sGenes in Human Adult Cornea scRNAseq

Please see the above directories for details of the methods and scripts.

# Datasets used

Multiple public datasets are used in the analysis including:

| Data          | Link          | Reference    |
| ------------- | ------------- | ------------- |
| CRF GWAS summary statistics  | [GWASCatalog Study GCST90308682](http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90308001-GCST90309000/GCST90308682/) | [1]
| CRF GWAS fine-mapping results using FINEMAP  | Supplementary Data 7 of the referred paper | [1]
| fine-mapping results using DAP-G for cis-e/sQTLs in all v8 GTEx tissues  | [https://zenodo.org/record/3517189#.Y-0UwcfP2Ul](https://zenodo.org/record/3517189#.Y-0UwcfP2Ul)  | [2]
| keratoconus GWAS summary statistics | Supplementary Data 15 of the referred paper | [3] |
| Human Adult Cornea scRNA-seq  | [http://retinalstemcellresearch.co.uk/CorneaCellAtlas](http://retinalstemcellresearch.co.uk/CorneaCellAtlas)  | [4] |

# Reference
[1] Jiang, X., Dellepiane, N., Pairo-Castineira, E. et al. Fine-mapping and cell-specific enrichment at corneal resistance factor loci prioritize candidate causal regulatory variants. *Commun Biol* 3, 762 (2020). [https://doi.org/10.1038/s42003-020-01497-w](https://doi.org/10.1038/s42003-020-01497-w) \
[2] Barbeira, A.N., Bonazzola, R., Gamazon, E.R. et al. Exploiting the GTEx resources to decipher the mechanisms at GWAS loci. *Genome Biol* 22, 49 (2021). [https://doi.org/10.1186/s13059-020-02252-4](https://doi.org/10.1186/s13059-020-02252-4) \
[3] Hardcastle, A.J., Liskova, P., Bykhovskaya, Y. et al. A multi-ethnic genome-wide association study implicates collagen matrix integrity and cell differentiation pathways in keratoconus. *Commun Biol* 4, 266 (2021). [https://doi.org/10.1038/s42003-021-01784-0](https://doi.org/10.1038/s42003-021-01784-0) \
[4] Collin, J., Queen, R., Zerti, D. et al. A single cell atlas of human cornea that defines its development, limbal progenitor cells and their interactions with the immune cells. *Ocul Surf* 21, 279 (2021). [https://doi: 10.1016/J.JTOS.2021.03.010](https://doi.org/10.1016/J.JTOS.2021.03.010)

# Citation
Jiang X, Boutin T and Vitart V (2023) Colocalization of corneal resistance factor GWAS loci with GTEx e/sQTLs highlights plausible candidate causal genes for keratoconus postnatal corneal stroma weakening. Front. Genet. 14:1171217. doi: 10.3389/fgene.2023.1171217
