library(tidyverse)
library(fda)
data <- read.csv("Data_Filter/data.csv")
G <- read.csv("Data_Filter/data_snp.csv", row.names = 1)
PTID <- row.names(G)
L <- table(data$PTID)[PTID]
G <- G[rep(PTID, L - 1), ]
data <- data %>% group_by(PTID) %>% mutate(diff = Y - first(Y))
data <- subset(data, select = -c(Y))
colnames(data)[which(names(data) == "diff")] <- "Y"
data.bl <- data %>% group_by(PTID) %>% slice(1)
data.bl <- data.frame(data.bl, row.names = 1)
X <- subset(data.bl, select = -c(Y))
data. <- data %>% group_by(PTID) %>% slice(2:n())
Y <- subset(data., select = c(PTID, AGE, Y))
PTID.Y <- Y$PTID
X <- X[PTID.Y, ]
colnames(X)[which(names(X) == "AGE")] <- "AGE_BL"
XY <- cbind(X, Y[, 2:3])
XY <- XY[rownames(G), ]
X <- subset(XY, select = -c(Y))
Y <- subset(XY, select = c(Y))
# .V <- list()
# for (i in 1:length(PTID)) {
#   X. <- X[grep(PTID[i], row.names(X)), , drop = FALSE]
#   elem <- X.$AGE - X.$AGE_BL
#   L. <- length(elem)
#   a <- matrix(rep(elem, L.), L., L.)
#   a[lower.tri(a)] <- t(a)[lower.tri(a)]
#   .V[[i]] <- solve(a)
# }
# V <- do.call(bdiag, .V)
# saveRDS(V, "Data_Filter/v.rds")
write.csv(X, "Data_Filter/x.csv", quote = FALSE)
write.csv(G, "Data_Filter/G.csv", quote = FALSE)
write.csv(Y, "Data_Filter/y.csv", quote = FALSE)