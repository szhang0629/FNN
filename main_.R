library(MASS)
library(fda)
library(Deriv)
source('Common/common.R')
source('Common/divide.R')
main <- function(seed, vari = 0, p = 200, HL = 1, fct = "polyno") {
  vari. <- c("index", "train", "test", "mse1", "mse2", "cor1", "cor2")
  if (vari == 0) {
    ipath. <- paste0("../4_Output/Real/")
    Y <- read.csv("../2_Data/Real/y.csv", as.is = TRUE)
    library(tidyverse)
    data <- as.data.frame(Y %>% group_by(PTID) %>% filter(n() > 1) %>%
                            mutate(loc = (AGE - 50) / 50) %>%
                            arrange(loc, .by_group = TRUE) %>%
                            mutate(loc0 = min(loc)) %>%
                            mutate(diff = Y - first(Y)) %>% slice(2:n()) %>%
                            ungroup() %>% mutate(Y = diff/sd(diff)))
    Y <- as.matrix(subset(data, select = Y))
    X <- as.matrix(subset(data, select = c(loc, loc0)))
    rownames(Y) <- data$PTID
    rownames(X) <- data$PTID
    
    G <- read.csv("../2_Data/Real/g.csv", row.names = 1)
    pos <- as.numeric(substring(colnames(G), 2))
    pos <- (pos - min(pos))/(max(pos) - min(pos)) * 0.9 + 0.05
    G <- as.matrix(G)
    
    G <- G[rownames(Y), ]
    lz <- 1
    loc <- NULL
    lambda. = 10^((2:6)/2)
    HU <- 10
  } else {
    ipath. <- paste0("../4_Output/", p, "/", vari, "/")
    n <- 250
    set.seed(seed)
    loc <- 0:19/20 + runif(20, 0.01, 0.04)
    X <- NULL
    source('Data/gdata.R')
    Data <- gData(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, p, loc, fct)
    list2env(Data, envir = environment())
    G <- G[rownames(Y), ]
    lz <- ncol(Y)
    lambda. <- 10^(-2:2)
    HU <- 40
  }
  if (!is.null(rownames(Y))) groups <- divide(rownames(Y), seed, "name")
  else groups <- divide(1:nrow(Y), seed)
  list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
  ipath <- paste0(ipath., "FLM1.csv")
  Output.(rbind(vari.), ipath)
  if (!(seed %in% (read.csv(ipath)$index))) {
    source("LM/flm.R")
    bb0 <- cbb(norder = 4, pos, 15)
    flm1.p <- flm1(Y.train, G.train, bb0, pos, X.train, lambda./10)
    Y.train. <- pred.flm1(G.train, bb0, pos, flm1.p$para, X.train)
    Y.test. <- pred.flm1(G.test, bb0, pos, flm1.p$para, X.test)
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, index = seed)[, vari.], digits = 4),
        ipath, pscn = TRUE)
  }
  if (is.numeric(loc) && length(loc) > 10) {
    ipath <- paste0(ipath., "FLM2.csv")
    Output.(rbind(vari.), ipath)
    if (!(seed %in% (read.csv(ipath)$index))) {
      source("LM/flm2.R")
      prt(format(cbind(Error.flm2(Y.train, G.train, Y.test, G.test, pos, loc),
                       index = seed)[, vari.], digits = 4), ipath, pscn = TRUE)
    }
  }
  ipath <- paste0(ipath., "NN", HL, ".csv")
  Output.(rbind(vari.), ipath)
  if (!(seed %in% (read.csv(ipath)$index))) {
    sapply(paste0("NN/", list.files("NN")), source)
    A <- c(rep(list(sigmoid), HL), list(linear))
    E.train <- cbind(X.train, G.train)
    E.test <- cbind(X.test, G.test)
    Bases <- c(list(ncol(E.train)), as.list(rep(HU, HL)), 
               list(if (!is.null(loc)) length(loc) else 1))
    nn.p <- NN.p(Y.train, E.train, Bases, A = A, lambda. = lambda.)
    Y.train. <- pred(nn.p$Ac, E.train, A)
    Y.test. <- pred(nn.p$Ac, E.test, A)
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, index = seed)[, vari.], digits = 4), ipath,
        pscn = TRUE)
  }
}