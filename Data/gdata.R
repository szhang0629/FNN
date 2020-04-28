gData <- function(index, vari = 0, n = 250, p = 500, fct. = "polyno", 
                  method = "N"){
  set.seed(index)
  if (vari > 0) {
    set.seed(index)
    loc <- 0:19/20 + runif(20, 0.01, 0.04)
    noise = switch(vari, 0.3, 0.6, 1.2)
    source('Data/data_input.R')
    source('Data/fun_output.R')
    source('Data/data_output.R')
    Data. <- Data.Input(n = n, p = p)
    list2env(Data., envir = environment())
    source(paste0("Data/", fct., ".R"))
    force(fct)
    f <- fct(Data.$G, Data.$pos, index, align = !is.matrix(loc))
    Y <- Data.Output(f, loc, noise)
    if (substr(method, 1, 1) == "F") {
      ipath <- paste0("../4_Output/", p, "/", vari, "/", method, ".csv")
      bb0 <- cbb(norder = 4, pos, 60)
      bb1 <- create.bspline.basis(norder = 4, nbasis = 40)
      bbz <- cbb(norder = 4, c(loc), 20)
      Bases <- c(list(bb0), rep(list(bb1), as.numeric(substring(
        method, nchar(method)))), list(bbz))
      Data <- list(G = Data.$G, X = NULL, Y = Y, pos = Data.$pos, loc = loc, 
                   lambda. = 10^(-3:1), Bases = Bases, ipath = ipath)
    } else 
      Data <- list(G = Data.$G, X = NULL, Y = Y, pos = Data.$pos, loc = loc, 
                   lambda. = 10^(-2:2), HU = 40, 
                   ipath. = paste0("../4_Output/", p, "/", vari, "/"))
  } else {
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
    
    G <- read.csv("../2_Data/Real/g.csv", row.names = 1)
    pos <- as.numeric(substring(colnames(G), 2))
    pos <- (pos - min(pos))/(max(pos) - min(pos)) * 0.9 + 0.05
    G <- as.matrix(G)
    if (substr(method, 1, 1) == "F") {
      ipath <- paste0("../4_Output/Real/", method, ".csv")
      bb0 <- cbb(norder = 4, pos, 15)
      bb1 <- create.bspline.basis(norder = 4, nbasis = 10)
      bbz <- cbb(norder = 4, c(Y[, 'loc'], Y[, 'loc0']), 10)
      Bases <- c(list(bb0), rep(list(bb1), as.numeric(substring(
        method, nchar(method)))), list(bbz))
      Data <- list(G = G, X = NULL, Y = Y, pos = pos, loc = NULL, 
                   lambda. = 10^(3:7), Bases = Bases, ipath = ipath)
    } else {
      X <- Y[, -1, drop = FALSE]
      Y <- Y[, 1, drop = FALSE]
      rownames(Y) <- data$PTID
      rownames(X) <- data$PTID
      G <- G[rownames(Y), ]
      Data <- list(G = G, X = X, Y = Y, pos = pos, loc = NULL, 
                   lambda. = 10^(0:4), HU = 10, ipath. = "../4_Output/Real/")
    }
  }
  return(Data)
}