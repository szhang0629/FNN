source("source.R")
real <- function(seed, D = 3, lambda. = 10^(-2:1)) {
  vari. <- c("method", "train", "test", "lambda", "j", "nl")
  ipath <- paste0("../4_Output/Real2/", seed, ".csv")
  Output.(rbind(vari.), ipath)
    
  Y <- read.csv("../2_Data/data_y.csv")
  colnames(Y)[colnames(Y) == "AGE"] <- "loc" 
  Y$loc <- 0.9 * Y$loc + 0.05
  Y$loc <- round(Y$loc, digits = 3)
  library(tidyverse)
  Y <- Y %>% group_by(PTID) %>% filter(n() > 1)
  Y <- Y %>% group_by(PTID) %>% arrange(loc, .by_group = TRUE)
  Y <- Y %>% group_by(PTID) %>% mutate(loc0 = min(loc))
  Y <- Y %>% group_by(PTID) %>% mutate(diff = Y - first(Y))
  Y <- subset(Y, select = -Y)
  colnames(Y)[which(names(Y) == "diff")] <- "Y"
  Y <- as.data.frame(Y %>% group_by(PTID) %>% slice(2:n()) %>% ungroup())
  Y$Y <- scale(exp(Y$Y))
  
  # groups <- divide(Y$PTID, seed, "name")
  groups <- list(train = unique(Y$PTID)[-seed], test = unique(Y$PTID)[seed])
  X <- as.matrix(read.csv("../2_Data/data_x.csv", row.names = 1))
  G <- as.matrix(read.csv("../2_Data/data_g.csv", row.names = 1))
  G[G != 0] <- 0
  pos <- read.csv("../2_Data/data_pos.csv")$POS
  list2env(split(mget(c('Y', 'X', 'G')), groups), envir = environment())
  
  A <- c(rep(list(sigmoid), D - 1), list(linear))
  method2 <- paste0("FN", D - 1)
  # if (!(method2 %in% (read.csv(ipath)$method))) {
  # la <- ceiling(ncol(G.train)*(1 - exp(-rankMatrix(G.train)/ncol(G.train))))
  la <- 15
  lz <- 15
  # l <- round(ncol(G.train)^(1/3) * 20^(2/3))
  l <- 15
  # if (!(method2 %in% (read.csv(ipath)$method))) {
  bb0 <- cbb(pos, nbasis = la)
  bb1 <- create.bspline.basis(norder = round(sqrt(l)), nbasis = l)
  
  bbz <- cbb(c(0, unlist(Y.train$loc), 1), nbasis = lz)
  Bases.fn <- c(list(bb0), rep(list(bb1), D - 1), list(bbz))
  error.fn <- Error.fnn(Y.train, X.train, G.train, Y.test, X.test,
                        G.test, pos, Bases.fn, A = A, lambda. = lambda.)
  # prt(cbind(error.fn, method = method2, nl = paste(la, l, lz))[, vari.], 
  #     ipath)
  # }
  output <- format(cbind(error.fn, method = method2, 
                         nl = paste(la, l, lz))[, vari.], digits = 4)
  write.table(output, ipath, append = T, quote = F, sep = ",", row.names = F, 
              col.names = F)
  print(read.csv(ipath), row.names = FALSE)
}
pen.fun <- function(bb, order = 2) {
  if (is.basis(bb)) {
    if (order == 2)
      # return(bb$nbasis*bsplinepen(bb, order)/1e7)
      # return(bsplinepen(bb, order)/(sum(diag(bsplinepen(bb, 0)))^2)/10^2)
      # return(bb$nbasis*bsplinepen(bb, order)/1e4)
      return(bsplinepen(bb, order)/1e2)
      # return(bb$nbasis*bsplinepen(bb, order) /
      #          sum(diag(bsplinepen(bb, 0))) /1e6)
    else 
      return(bb$nbasis*bsplinepen(bb, order))
      # return(bb$nbasis*bsplinepen(bb, order)/sum(diag(bsplinepen(bb, 0))))
  } else return(1)
}
f2m. <- function(bb, loc) {
  if (!is.data.frame(loc)) {
    coef <- bb$nbasis^(1/2)
    # coef <- 1/sum(diag(bsplinepen(bb, 0)))
    # coef <- (bb$nbasis/sum(diag(bsplinepen(bb, 0))))^(1/2)
    B <- eval.basis(loc, bb) * coef
  } else if (length(loc) == 1 && is.data.frame(loc))
    B <- f2m.(bb, unlist(loc))
  else B <- (f2m.(bb, loc$loc) - f2m.(bb, loc$loc0))
  return(B)
}