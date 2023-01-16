library(EWCE)
library(ggplot2)
set.seed(1)

# create CellTypeDataset
mat <- read.table("../../data/scrnaseq/exprMatrix.tsv.gz", sep = '\t',header=TRUE)
genes <- mat[, 1]
genes <- gsub(".+[|]", "", genes)
mat <- data.frame(mat[, -1], row.names = genes)
mat <- as.matrix(mat)

## normalise expression data
mat <- EWCE::sct_normalize(mat)
## read celltypes
meta <- read.table(
  "../../data/scrnaseq/meta.tsv",
  header = T,
  sep = "\t",
  as.is = T,
  row.names = 1)

## create CellTypeDataset
ctd_path <- generate_celltype_data(
  exp = mat,
  annotLevels = list(meta$cluster),
  input_species = "human",
  output_species = "human",
  groupName = 'AdultCornea',
  savePath = tempdir())

ctd <- EWCE::load_rdata(ctd_path)

# write specificity, mean exp, and specificity quantiles
write.table(as.data.frame(as.matrix(ctd[[1]]$specificity)), file='ewce_specificity.tsv',sep='\t',quote = FALSE, row.names = FALSE)
write.table(as.data.frame(as.matrix(ctd[[1]]$mean_exp)), file='ewce_mean_exp.tsv',sep='\t',quote = FALSE, row.names = FALSE)
write.table(as.data.frame(as.matrix(ctd[[1]]$specificity_quantiles)), file='ewce_specificity_quantiles.tsv',sep='\t',quote = FALSE, row.names = FALSE)

# gene list
esgenes <- read.table("./input/esgenes.txt")
esgenes <- esgenes$V1
nearest_genes <- read.table("./input/nearest_genes.txt")
nearest_genes <- nearest_genes$V1
nearest_genes_exclude_esgenes <- read.table("./input/nearest_genes_exclude_esgenes.txt")
nearest_genes_exclude_esgenes <- nearest_genes_exclude_esgenes$V1
esgenes_nearest_substitute <- read.table("./input/esgenes_nearest_substitute.txt")
esgenes_nearest_substitute <- esgenes_nearest_substitute$V1

# run enrichment tests
reps <- 10000
esgenes_enrich <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = esgenes, 
                                                reps = reps,
                                                annotLevel = 1)


nearest_genes_enrich <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                  sctSpecies = "human",
                                                  genelistSpecies = "human",
                                                  hits = nearest_genes, 
                                                  reps = reps,
                                                  annotLevel = 1)

nearest_genes_exclude_esgenes_enrich <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                  sctSpecies = "human",
                                                  genelistSpecies = "human",
                                                  hits = nearest_genes_exclude_esgenes,
                                                  reps = reps,
                                                  annotLevel = 1)
                                          
esgenes_nearest_substitute_enrich <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                  sctSpecies = "human",
                                                  genelistSpecies = "human",
                                                  hits = esgenes_nearest_substitute,
                                                  reps = reps,
                                                  annotLevel = 1)

write.table(nearest_genes_exclude_esgenes_enrich$results, file='nearest_genes_exclude_esgenes_enrich.tsv',sep='\t',quote = FALSE, row.names = FALSE)
write.table(esgenes_enrich$results, file='esgenes_enrich.tsv',sep='\t',quote = FALSE, row.names = FALSE)
write.table(nearest_genes_enrich$results, file='nearest_genes_enrich.tsv',sep='\t',quote = FALSE, row.names = FALSE)
write.table(esgenes_nearest_substitute_enrich$results, file='esgenes_nearest_substitute_enrich.tsv',sep='\t',quote = FALSE, row.names = FALSE)


# non-conditional
non_conditionalres <- rbind(
  data.frame(esgenes_enrich$results,
             list="e/sgenes"),
  data.frame(nearest_genes_enrich$results,
             list="nearest genes"),
  data.frame(nearest_genes_exclude_esgenes_enrich$results,
             list="nearest genes excluding e/sgenes"),
  data.frame(esgenes_nearest_substitute_enrich$results,
             list = "e/sgenes & nearest genes")) #i.e. esgenes nearest genes substitute

figure1 <- EWCE::ewce_plot(total_res = non_conditionalres, mtc_method = "BH")
figure1 <- figure1$plain
figure1 <- figure1 + theme(axis.text.x = element_text(angle = 70))+coord_cartesian(ylim = c(0,10),expand = FALSE)+ scale_y_continuous(breaks = seq(0, 10, by = 3))


# conditional analysis
# esgenes
esgenes_enrich_cond_cssc <- EWCE:: bootstrap_enrichment_test(
  sct_data = ctd,
  hits = esgenes,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = reps,
  annotLevel = 1,
  controlledCT = "corneal_stromal_stem_cells")

esgenes_enrich_cond_sck <- EWCE:: bootstrap_enrichment_test(
  sct_data = ctd,
  hits = esgenes,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = reps,
  annotLevel = 1,
  controlledCT = "corneal_stroma_keratocytes")

#nearest_genes
nearest_genes_cond_cssc <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                                   sctSpecies = "human",
                                                                   genelistSpecies = "human",
                                                                   hits = nearest_genes, 
                                                                   reps = reps,
                                                                   annotLevel = 1,
                                                                   controlledCT = "corneal_stromal_stem_cells")
nearest_genes_cond_csk <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                                  sctSpecies = "human",
                                                                  genelistSpecies = "human",
                                                                  hits = nearest_genes, 
                                                                  reps = reps,
                                                                  annotLevel = 1,
                                                                  controlledCT = "corneal_stroma_keratocytes")
nearest_genes_cond_lf <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                                 sctSpecies = "human",
                                                                 genelistSpecies = "human",
                                                                 hits = nearest_genes, 
                                                                 reps = reps,
                                                                 annotLevel = 1,
                                                                 controlledCT = "limbal_fibroblasts")
nearest_genes_cond_lsk <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                                  sctSpecies = "human",
                                                                  genelistSpecies = "human",
                                                                  hits = nearest_genes, 
                                                                  reps = reps,
                                                                  annotLevel = 1,
                                                                  controlledCT = "limbal_stroma_keratocytes")

#esgenes_nearest_substitute
esn_sub_cond_cssc <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                           sctSpecies = "human",
                                                           genelistSpecies = "human",
                                                           hits = esgenes_nearest_substitute, 
                                                           reps = reps,
                                                           annotLevel = 1,
                                                           controlledCT = "corneal_stromal_stem_cells")
esn_sub_cond_csk <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                          sctSpecies = "human",
                                                          genelistSpecies = "human",
                                                          hits = esgenes_nearest_substitute, 
                                                          reps = reps,
                                                          annotLevel = 1,
                                                          controlledCT = "corneal_stroma_keratocytes")
esn_sub_cond_lf <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                         sctSpecies = "human",
                                                         genelistSpecies = "human",
                                                         hits = esgenes_nearest_substitute, 
                                                         reps = reps,
                                                         annotLevel = 1,
                                                         controlledCT = "limbal_fibroblasts")
esn_sub_cond_lsk <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                          sctSpecies = "human",
                                                          genelistSpecies = "human",
                                                          hits = esgenes_nearest_substitute, 
                                                          reps = reps,
                                                          annotLevel = 1,
                                                          controlledCT = "limbal_stroma_keratocytes")


all_conditionalres <- rbind(
  data.frame(esgenes_enrich_cond_cssc$results,
             list="e/sgenes, corneal stromal stem cells controlled"),
  data.frame(esgenes_enrich_cond_sck$results,
             list="e/sgenes, corneal stroma keratocytes controlled"),
  data.frame(nearest_genes_cond_cssc$results,
             list="nearest genes, corneal stromal stem cells controlled"),
  data.frame(nearest_genes_cond_csk$results,
             list="nearest genes, corneal stroma keratocytes controlled"),
  data.frame(nearest_genes_cond_lf$results,
             list="nearest genes, limbal fibroblasts controlled"),
  data.frame(nearest_genes_cond_lsk$results,
             list="nearest genes, limbal stroma keratocytes controlled"),
  data.frame(esn_sub_cond_cssc$results,
             list="e/sgenes & nearest genes, corneal stromal stem cells controlled"),
  data.frame(esn_sub_cond_csk$results,
             list="e/sgenes & nearest genes, corneal stroma keratocytes controlled"),
  data.frame(esn_sub_cond_lf$results,
             list="e/sgenes & nearest genes, limbal fibroblasts controlled"),
  data.frame(esn_sub_cond_lsk$results,
             list="e/sgenes & nearest genes, limbal stroma keratocytes controlled")) #i.e. esgenes nearest genes substitute

figure <- EWCE::ewce_plot(total_res = all_conditionalres, mtc_method = "BH")
figure <- figure$plain
figure <- figure + theme(axis.text.x = element_text(angle = 70))+coord_cartesian(ylim = c(0,10),expand = FALSE)+ scale_y_continuous(breaks = seq(0, 10, by = 3))

#save image
save.image("ewce.RData")
