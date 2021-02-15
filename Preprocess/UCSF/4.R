data <- read.csv("Data_Filter/data.csv")
fct.inp <- read.csv("Data_Filter/data_snp.csv", row.names = 1)
X <- data[, c("PTID", "PTGENDER", "PTEDUCAT")]
X <- unique(X)
rownames(X) <- X$PTID
sca.inp <- X[, c("PTGENDER", "PTEDUCAT")]
fct.out <- data[, c("PTID", "AGE", "Y")]
fct.out$Y <- (fct.out$Y - mean(fct.out$Y)) / sd(fct.out$Y)



write.csv(fct.inp, "Data_Filter/fct_inp.csv", quote = FALSE)
write.csv(sca.inp, "Data_Filter/sca_inp.csv", quote = FALSE)
write.csv(Y, "Data_Filter/fct_out.csv", quote = FALSE)
