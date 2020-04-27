## ----------Data Substraction----------
Data.Input <- function(n, p){
  # ----------Substract data given N and p----------
  maf.interval <- c(0.01,0.5)
  variant <- function(x){length(unique(x)) > 1}
  ### Get data for use
  ped <- read.table("../2_Data/ped.ped", sep = "\t")
  info <- read.table("../2_Data/info.info", header = T)
  # ped: 1092 obs. of 12735 variables
  # info: 12735 obs. of 4 variables
  # n = 1092, p = 12735 high dimensional question
  
  ## N samples chosen among 1092 objects for simulation
  smp.idx <- sample(1:nrow(ped), n) # smp as sample order
  ## maf interval SNP index
  maf.idx <- (info$maf > maf.interval[1] & info$maf < maf.interval[2])
  geno <- ped[smp.idx, maf.idx]
  pos <- info$pos
  loc <- pos[maf.idx]
  
  ## Delete void data
  vrt <- apply(geno,2,variant)  # see variability with a genomic region
  loc <- loc[vrt]
  geno <- as.matrix(geno[,vrt]); # get rid of individuals with no variability
  
  ## Truncated SNP index
  seg.pos <- sample(1:(length(loc) - p + 1), 1)
  idx.trun <- seg.pos:(seg.pos + p - 1)
  geno <- geno[, idx.trun]
  loc <- loc[idx.trun]
  
  G <- as.matrix(geno) # get rid of individuals with no variability
  # rownames(G) <- 1:nrow(G)
  pos <- (loc - loc[1])/(loc[length(loc)] - loc[1])
  
  return(list(G = G, pos = pos))
}
