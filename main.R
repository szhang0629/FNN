source("source.R")
main <- function(seed, vari = 0, p = 200, D = 2, fct = "polyno") {
  vari. <- c("index", "train", "test", "mse1", "mse2", "cor1", "cor2")
  if (vari == 0) {
    ipath <- paste0("../4_Output/Real2/FN", (D - 1), ".csv")
    Y <- read.csv("../2_Data/data_y.csv", as.is = TRUE)
    colnames(Y)[colnames(Y) == "AGE"] <- "loc"
    
    bbz <- create.bspline.basis(norder = 4, nbasis = 15)
    Y$loc <- loc.scale(Y$loc, bbz)
    
    library(tidyverse)
    Y <- Y %>% group_by(PTID) %>% filter(n() > 1) %>% 
      arrange(loc, .by_group = TRUE) %>% mutate(loc0 = min(loc)) %>% 
      mutate(diff = Y - first(Y)) %>% slice(2:n()) %>% ungroup()
    Y <- as.data.frame(Y) 
    Y <- subset(Y, select = -Y)
    colnames(Y)[which(names(Y) == "diff")] <- "Y"
    Y$Y <- Y$Y/sd(Y$Y)
    # X <- read.csv("../2_Data/data_x_APOE.csv", row.names = 1)
    # X <- subset(X, select = -APOE4)
    # X <- as.matrix(cbind(subset(Y, select = -c(PTID, Y)), X[Y$PTID, ]))
    X <- as.matrix(read.csv("../2_Data/data_x_APOE.csv", row.names = 1))
    # X <- subset(X, select = -APOE4)
    G <- as.matrix(read.csv("../2_Data/data_g.csv", row.names = 1))
    pos <- read.csv("../2_Data/data_pos.csv")$POS
    loc <- NULL
    lambda. = 10^(-1:3)
  } else {
    ipath <- paste0("../4_Output/", p, "/", vari, "/FN", D - 1, ".csv")
    n <- 250
    set.seed(seed)
    # loc <- data.frame(PTID = rep(1:n, each = 20), loc = runif(20*n))
    # loc <- data.frame(PTID = rep(1:n, each = 20), loc = rep(runif(20), n))
    loc <- runif(20)
    Data <- gData(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, p, loc, fct)
    list2env(Data, envir = environment())
    rownames(G) <- 1:n
    if (is.data.frame(loc)) {
      bbz <- cbb(norder = 4, Y$loc, 20)
      # Y$loc <- loc.scale(Y$loc, bbz)
    } else {
      bbz <- cbb(norder = 4, loc, 20)
      # loc <- loc.scale(loc, bbz)
    }
    X <- NULL
    lambda. = list(1, 10, 100, 1000)
  }
  if (is.data.frame(Y)) groups <- divide(Y$PTID, seed, "name")
  else groups <- divide(1:nrow(Y), seed)
  list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
  Output.(rbind(vari.), ipath)
  # if (!(seed %in% (read.csv(ipath)$index))) {
    A <- c(rep(list(sigmoid), D - 1), list(linear))
    # la <- ceiling(ncol(G.train)*(1 - exp(-rankMatrix(G.train)/ncol(G.train))))
    la <- 20
    l <- ceiling(sqrt(bbz$nbasis * la))
    # bb0 <- create.bspline.basis(norder = 4, nbasis = la)
    bb0 <- cbb(norder = 4, pos, la)
    bb1 <- create.bspline.basis(norder = 4, nbasis = l)
    Bases <- c(list(bb0), rep(list(bb1), D - 1), list(bbz))
    
    fnn.p <- FNN.p(Y.train, X.train, G.train, Bases, A, lambda., pos, loc)
    if (is.data.frame(Y.train)) {
      Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos, 
                       loc = subset(Y.train, select = -Y))
      Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos, 
                      loc = subset(Y.test, select = -Y))
    } else {
      Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos, loc)
      Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos, loc)
    }
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, index = seed)[, vari.], digits = 4), ipath, pscn = TRUE)
  # }
}