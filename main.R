library(fda)
library(Deriv)
library(MASS)
source('Common/common.R')
source('Common/divide.R')
main <- function(seed, vari = 0, HL = 1, p = 200, fct = "polyno") {
  vari. <- c("index", "train", "test", "mse1", "mse2", "cor1", "cor2")
  if (vari == 0) {
    ipath <- paste0("../4_Output/Real/FN", HL, ".csv")
    Y <- read.csv("../2_Data/Real/y.csv", as.is = TRUE)
    # Y$loc <- loc.scale(Y$loc, bbz)
    library(tidyverse)
    data <- as.data.frame(Y %>% group_by(PTID) %>% filter(n() > 1) %>%
                            mutate(loc = (AGE - 50) / 50) %>%
                            arrange(loc, .by_group = TRUE) %>%
                            mutate(loc0 = min(loc)) %>%
                            mutate(diff = Y - first(Y)) %>% slice(2:n()) %>%
                            ungroup() %>% mutate(Y = diff/sd(diff)))
    Y <- as.matrix(subset(data, select = c(Y, loc, loc0)))
    rownames(Y) <- data$PTID
    X <- NULL
    
    G <- read.csv("../2_Data/Real/g.csv", row.names = 1)
    pos <- as.numeric(substring(colnames(G), 2))
    pos <- (pos - min(pos))/(max(pos) - min(pos)) * 0.9 + 0.05
    G <- as.matrix(G)
    
    loc <- NULL
    lambda. = 10^((6:10)/2)
    # lambda. = list(c(10^3, 10^3, 10^4))
    bb0 <- cbb(norder = 4, pos, 15)
    bb1 <- create.bspline.basis(norder = 4, nbasis = 10)
    bbz <- cbb(norder = 4, c(Y[, 'loc'], Y[, 'loc0']), 10)
    Bases <- c(list(bb0), rep(list(bb1), HL), list(bbz))
  } else {
    ipath <- paste0("../4_Output/", p, "/", vari, "/FN", HL, ".csv")
    n <- 250
    set.seed(seed)
    loc <- 0:19/20 + runif(20, 0.01, 0.04)
    source('Data/gdata.R')
    Data <- gData(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, p, loc, fct)
    list2env(Data, envir = environment())
    # Y <- cbind(Y, loc)
    # loc <- NULL
    X <- NULL
    lambda. = 10^(-3:1)
    bb0 <- cbb(norder = 4, pos, 60)
    bb1 <- create.bspline.basis(norder = 4, nbasis = 40)
    bbz <- cbb(norder = 4, c(loc), 20)
    Bases <- c(list(bb0), rep(list(bb1), HL), list(bbz))
  }
  if (!is.null(rownames(Y))) groups <- divide(rownames(Y), seed, "name")
  else groups <- divide(1:nrow(Y), seed)
  list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
  Output.(rbind(vari.), ipath)
  # if (!(seed %in% (read.csv(ipath)$index))) {
    sapply(paste0("FNN/", list.files("FNN")), source)
    A <- c(rep(list(sigmoid), HL), list(linear))
    fnn.p <- FNN.p(Y.train, X.train, G.train, Bases, A, lambda., pos, loc)
    if (is.null(loc)) {
      Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos, 
                       Y.train[, -1, drop = FALSE])
      Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos, 
                      Y.test[, -1, drop = FALSE])
      error <- Error(Y.train[, 1, drop = FALSE], Y.test[, 1, drop = FALSE], 
                     Y.train., Y.test.)
    } else {
      Y.train. <- pred(fnn.p$Ac, X.train, G.train, A, Bases, pos, loc)
      Y.test. <- pred(fnn.p$Ac, X.test, G.test, A, Bases, pos, loc)
      error <- Error(Y.train, Y.test, Y.train., Y.test.)
    }
    prt(format(cbind(error, index = seed)[, vari.], digits = 4), ipath, 
        pscn = TRUE)
  # }
}