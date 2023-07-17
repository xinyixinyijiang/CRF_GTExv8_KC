# scrnaseq_celltype_enrichment

`./ewce.R` is the main script for running enrichment of CRF signals linked GTEx v8 e/sGenes in [Human Adult Cornea scRNAseq](https://doi.org/10.1016/J.JTOS.2021.03.010). `./input/` contains the input gene lists for `./ewce.R`.

## Environment: 
- [EWCE](https://nathanskene.github.io/EWCE/articles/EWCE.html) docker
    - R: 4.2.1
    - image: neurogenomicslab/ewce:1.5.7
    - DIGEST: 3bbc3293c95c
    - [EWCE DOCKER HUB](https://hub.docker.com/r/neurogenomicslab/ewce)