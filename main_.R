source('source_.R')
main <- function(seed, vari = 0, p = 200, D = 3, lambda. = 10^(-2:2)) {
  vari. <- c("index", "train", "test", "cor1", "cor2", "lambda", "j")
  if (vari == 0) {
    ipath. <- paste0("../4_Output/Real2/")
    Y <- read.csv("../2_Data/data_y.csv", as.is = TRUE)
    colnames(Y)[colnames(Y) == "AGE"] <- "loc"
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
    X <- as.matrix(subset(Y, select = -c(PTID, Y)))
    rownames(X) <- Y$PTID
    G <- as.matrix(read.csv("../2_Data/data_g.csv", row.names = 1))
    G <- G[Y$PTID, ]
    pos <- read.csv("../2_Data/data_pos.csv")$POS
    Y <- subset(Y, select = c(PTID, Y))
    lz <- 1
  } else {
    ipath. <- paste0("../4_Output/", p, "/", vari, "/")
    n <- 250
    set.seed(seed)
    # loc <- data.frame(PTID = rep(1:n, each = 20), loc = rep(runif(1:20), n))
    loc <- runif(20)
    Data <- gData(seed, noise = switch(vari, 0.3, 0.6, 1.2), n, p, loc, "linear")
    list2env(Data, envir = environment())
    rownames(G) <- 1:n
    # G <- G[Y$PTID, ]
    X <- NULL
    # Y <- subset(Y, select = c(PTID, Y))
    lz <- ncol(Y)
  }
  if (is.data.frame(Y)) groups <- divide(Y$PTID, seed, "name")
  else groups <- divide(1:nrow(Y), seed)
  list2env(sep(mget(c('Y', 'X', 'G')), groups), envir = environment())
  ipath <- paste0(ipath., "FLM1.csv")
  Output.(rbind(vari.), ipath)
  if (!(seed %in% (read.csv(ipath)$index))) {
    bb0 <- create.bspline.basis(norder = 4, nbasis = 50)
    flm1.p <- flm1(Y.train, G.train, Y.test, G.test, bb0, pos, X.train, 
                   X.test, lambda. = c(10^(-3:1)))
    Y.train. <- pred.flm1(G.train, bb0, pos, flm1.p$para, X.train)
    Y.test. <- pred.flm1(G.test, bb0, pos, flm1.p$para, X.test)
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, index = seed, j = 0, 
                     lambda = flm1.p$lambda)[, vari.], digits = 4), ipath)
    # print(read.csv(ipath), row.names = FALSE)
  }
  ipath <- paste0(ipath., "NN", D - 1, ".csv")
  Output.(rbind(vari.), ipath)
  # if (!(seed %in% (read.csv(ipath)$index))) {
    A <- c(rep(list(sigmoid), D - 1), list(linear))
    E.train <- cbind(X.train, G.train)
    E.test <- cbind(X.test, G.test)
    Bases <- c(list(ncol(E.train)), as.list(rep(15, D - 1)), list(lz))
    nn.p <- NN.p(Y.train, E.train, Bases, A = A, lambda. = lambda.)
    Y.train. <- pred(nn.p$Ac, E.train, A)
    Y.test. <- pred(nn.p$Ac, E.test, A)
    error <- Error(Y.train, Y.test, Y.train., Y.test.)
    prt(format(cbind(error, index = seed, lambda = nn.p$lambda, 
                     j = nn.p$j)[, vari.], digits = 4), ipath)
    # print(read.csv(ipath), row.names = FALSE)
  # }
}