library(tidyverse)
data <- read.csv("Data_Filter/data.csv")
data.g <- read.csv("Data_Filter/data_snp.csv", row.names = 1)
POS <- read.csv("Data_Filter/pos.csv")
## Transform Data
pos <- (POS - min(POS))/(max(POS) - min(POS))
## Seperate X
# X <- data[, c("PTID", "PTGENDER", "PTEDUCAT", "APOE4")]
X <- data[, c("PTID", "PTGENDER", "PTEDUCAT")]
X <- data.frame(unique(X))
X$PTEDUCAT <- scale(X$PTEDUCAT)
rownames(X) <- X$PTID
X <- as.matrix(X[, -1])
## Seperate Y and loc
data.y <- data[, c("PTID", "AGE", "Y")]
data.y <- data.frame(data.y)
PTID <- rownames(X)
data.g <- data.g[PTID, ]
loc <- list()
Y <- list()
for (i in 1:length(PTID)) {
  loc. <- data.y[data.y$PTID == PTID[i], "AGE"]
  loc[[i]] <- loc.
  Y. <- data.y[data.y$PTID == PTID[i], "Y"]
  Y[[i]] <- Y.
}
names(loc) <- PTID
names(Y) <- PTID
saveRDS(list(G = data.g, pos = pos, X = X, Y = Y, loc = loc),"Data_Func/R.rds")
