readVCF <- function(vcf) {
  library(vcfR)
  library(Matrix)
  
  ## load VCF
  vc <- read.vcfR(vcf)
  if (nrow(vc) == 0L)
    stop('null genotype ', vcf)
  gn <- sub('.vcf.gz', '', basename(vcf))
  
  ## extract non-synonymous variants
  gt <- extract.gt(vc)
  
  ## sample names, counts
  sbj <- colnames(gt)
  
  ## genomic map
  mp <- vcfR2tidy(vc, T, info_fields = '.')$fix[, 1:5]
  names(mp)[1] <- 'CHR'
  ## chromosome map
  cm <- 1:26
  names(cm) <- c(1:22, 'X', 'Y', 'M', 'XY')
  mp <- within(mp, CHR <- cm[CHR])
  
  ## variant IDs
  gvr <- sprintf('V%05X', 1L:nrow(mp))
  
  ## the GT to dosage dictionary.
  dc <- c(
    0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
    1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
    1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
    1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
    1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)
  ## bi-allelic format to dosage format
  gt <- as.integer(gsub('[/|]', '', gt))
  gt <- dc[gt + 1L]
  gt <- Matrix(gt, length(gvr), length(sbj), dimnames = list(gvr = gvr,
                                                           sbj = sbj))
  
  rt <- within(list(),
               {
                 map <- mp
                 gmx <- gt
                 bp1 <- tail(mp$POS, 1L)
                 bp0 <- head(mp$POS, 1L)
                 chr <- mp$CHR[1]
                 vcf <- vcf
               })
  
  invisible(rt)
}
