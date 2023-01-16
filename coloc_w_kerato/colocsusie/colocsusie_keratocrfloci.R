library(coloc)
library(stringr)
library(dplyr)

crfloci <- read.table('../loci/CRF115loci.tsv', sep='\t', header=TRUE)
colnames(crfloci) <- c('loci', 'chr','start','end','starthg38', 'endhg38')
kerato_crfloci <- read.table('../crf_kerato_overlaploci.tsv', sep='\t', header=FALSE)
colnames(kerato_crfloci) <- c(
  'crfloci_chr','crfloci_start','crfloci_end','crfloci', 
  'keratoloci_chr', 'keratoloci_start', 'keratoloci_end', 'keratoloci')

#### FUNCTIONS #####
prepare_datasets <- function(datapath, ldpath, sdy, nsample, type){
  data <- read.table(datapath, sep = '\t', header = TRUE)
  R <- as.matrix(read.csv(ldpath, sep = ' ', header = FALSE))
  R <- R[,colSums(is.na(R))<nrow(R)]
  dimnames(R) <- list(data$rsid, data$rsid)
  if (is.null(sdy)){
    colocds <- list(data$beta,data$se * data$se, nsample,type, data$maf, R, data$rsid, data$position)
    names(colocds) <- c("beta","varbeta", "N","type", "MAF", "LD", "snp", "position")
  }else{
    colocds <- list(data$beta,data$se*data$se, nsample,sdy, type, data$maf, R, data$rsid, data$position)
    names(colocds) <- c("beta","varbeta", "N","sdY","type", "MAF", "LD", "snp", "position")
  }
  if (is.null(check_dataset(colocds, req='LD'))){
    return(colocds) 
  }else{
    return('check_dataset not equal to NULL')
  }
}

output_susieres <- function(datapath, susieres, outputpath, outputname){
  data <- read.table(datapath, sep = '\t', header = TRUE)
  pips <- data.frame(summary(susieres)[[1]])
  pips <- pips[pips$cs > 0, ]
  pips <- merge(pips, data, by.x = 0, by.y = 0)
  pips$Row.names <- NULL
  pips$variable <- NULL
  colnames(pips) <- c("pip", "cs", "rsid", "position", "beta", "se", "maf")
  pips <- pips[order(pips$pip, decreasing = TRUE), ]
  write.table(
    pips,
    file = paste0(outputpath, outputname, ".susieoutput"),
    quote = FALSE, sep = "\t", row.names = FALSE)
}

#### CRF KERATO COLOC #####
colocres <- NULL
for (i in as.numeric(row.names(kerato_crfloci))){
  chrom <- kerato_crfloci$crfloci_chr[[i]]
  start <- kerato_crfloci$crfloci_start[[i]]
  end <- kerato_crfloci$crfloci_end[[i]]
  crflociid <- str_split(kerato_crfloci$crfloci[[i]], 'locus')[[1]][2]
  keratolociid <- str_split(kerato_crfloci$keratoloci[[i]], 'locus')[[1]][2]
  ds1 <- c(
    paste0("../input/crf/crf",crflociid,".colocinput"),
    paste0("../input/crf/crf",crflociid,".ld"),
    "../output/susie/crf/",
    paste0("crf",crflociid)
  )
  
  ds2 <- c(
    paste0("../input/kerato/kerato", keratolociid, "_crf", crflociid,".colocinput"),
    paste0("../input/kerato/kerato", keratolociid, "_crf", crflociid,".ld"),
    "../output/susie/kerato/",
    paste0("kerato",keratolociid,"_crf", crflociid)
  )
  
  crfsdy <- 1.8233759482860157
  crf <- prepare_datasets(ds1[1], ds1[2], crfsdy, 72301, "quant")
  kerato <- prepare_datasets(ds2[1], ds2[2], sdy=NULL, 121216, "cc")
  scrf <- runsusie(crf)
  skerato <- runsusie(kerato)
  output_susieres(ds1[1], scrf, ds1[3], ds1[4])
  output_susieres(ds2[1], skerato, ds2[3], ds2[4])
  if (length(scrf$sets$cs) != 0 && length(skerato$sets$cs) != 0) {
    coloc.res <- coloc.susie(scrf,skerato)
    res <- as.data.frame(coloc.res$summary)
    res$crfloci <- crflociid
    res$keratoloci <- keratolociid
    colocres <- bind_rows(colocres, res)
    print(paste(as.character(i), 'done'))
  }else{
    if (length(scrf$sets$cs) != 0){
      print(paste('overlapping loci, kerato', as.character(keratolociid), 'failed', '(crf loci',crflociid,'), because no cs found by susie in keratoconus GWAS'))
    }else{
      print(paste('overlapping loci, kerato', as.character(keratolociid), 'failed', '(crf loci',crflociid,'), because no cs found by susie in crf GWAS'))
    }
  }
}

write.table(colocres, '../output/coloc/crf_kerato/colocres.tsv', sep='\t', quote=FALSE, row.names=FALSE)
