library(tidyverse)
library(lubridate)
source('readVCF.R')
## Phenotype of LR
data.phe <- read.csv('Data_Original/UCSFFSX51_08_01_16.csv', as.is = TRUE)
# data.phe <- data.phe[, c("RID", "EXAMDATE", 'VERSION', "LHIPQC", "ST29SV")]
data.phe <- data.phe[, c("RID", "EXAMDATE", 'VERSION', "RHIPQC", "ST88SV")]
data.phe <- subset(data.phe[data.phe[, "RHIPQC"] == "Pass", ], 
                   select = -c(RHIPQC))
colnames(data.phe)[ncol(data.phe)] <- "Y"

data.phe$EXAMDATE <- as.Date(data.phe$EXAMDATE, format = c("%m/%d/%Y"))
data.phe$VERSION <- as.Date(data.phe$VERSION, format = c("%m/%d/%Y"))
data.phe <- data.phe %>% group_by(RID, EXAMDATE) %>% 
  filter(VERSION == max(VERSION)) %>% ungroup()
data.phe <- subset(data.phe, select = -c(VERSION))

data.phe <- data.phe %>% group_by(RID, EXAMDATE) %>% slice(n())

# data.phe <-  data.phe %>% group_by(RID) %>% filter(n() > 1)
## Covariate
data.aux <- read.csv('Data_Original/ADNIMERGE.csv', as.is = TRUE)
data.aux <- data.aux[data.aux$VISCODE == 'bl', ]
data.aux <- data.aux[, c('RID', 'PTID', 'EXAMDATE', 'AGE',
                         'PTGENDER', 'PTEDUCAT', 'APOE4')]
# data.aux <- data.aux[, c('RID', 'PTID', 'EXAMDATE', 'AGE', 'PTGENDER', 
#                          'PTEDUCAT')]
data.aux$EXAMDATE <- as.Date(data.aux$EXAMDATE)
colnames(data.aux)[colnames(data.aux) == "EXAMDATE"] <- 'EXAMDATE_bl'
colnames(data.aux)[colnames(data.aux) == "AGE"] <- 'AGE_bl'
data.aux$PTEDUCAT <- scale(data.aux$PTEDUCAT)
## Combine Phenotype and Covariate
data <- merge(data.aux, data.phe, by = 'RID')
data$AGE <- time_length(difftime(data$EXAMDATE, data$EXAMDATE_bl), "years") + 
  data$AGE_bl
data <- subset(data, select = -c(EXAMDATE_bl, AGE_bl, EXAMDATE))
data <- data[order(data$RID, data$AGE), ]
levels(data$PTGENDER) <- c(0, 1)
data$PTGENDER[data$PTGENDER == "Male"] <- 1
data$PTGENDER[data$PTGENDER == "Female"] <- 0
# age <- data$AGE
# data$AGE <- (data$AGE - min(age))/(max(age) - min(age))
# data$Y <- scale(log(data$Y))
# data$Y <- (data$Y)/100
# data$Y <- c(scale(data$Y))
data <- subset(data, select = -c(RID))
## Genotype
data.vcf <- readVCF("Data_Original/APOE.vcf.gz")
data.gmx <- data.vcf$gmx
data.snp <- t(as.matrix(data.gmx))
colnames(data.snp) <- data.vcf$map$ID
data.POS <- data.vcf$map$POS
## Select Common Samples
PTID. <- intersect(data$PTID, rownames(data.snp))
data.snp. <- data.snp[PTID., ]
data <- data[data$PTID %in% PTID., ]
## Remove columns of 10 percent NAs, Otherwise, use colMeans to substitute
na.col <- colMeans(is.na(data.snp.))
data.g <- data.snp.[, na.col < 0.1]
POS <- data.POS[na.col < 0.1]
for (i in 1:ncol(data.g))
  data.g[is.na(data.g[,i]), i] <- mean(data.g[,i], na.rm = TRUE)
variant <- function(x){length(unique(x)) > 1}
vrt <- apply(data.g, 2, variant)
POS <- POS[vrt]
data.g <- as.matrix(data.g[ ,vrt])
data.y <- data[, c('PTID', 'AGE', 'Y')]
data.x <- unique(data[, c('PTID', 'PTGENDER', 'PTEDUCAT', 'APOE4')])
# data.x <- unique(data[, c('PTID', 'PTGENDER', 'PTEDUCAT')])
data.x <- data.x %>% remove_rownames %>% column_to_rownames(var = "PTID")
POS <- (POS - min(POS)) / (max(POS) - min(POS))
write.csv(data.x, "Data_N/data_x_APOE.csv", quote = FALSE)
# write.csv(data.x, "Data_N/data_x.csv", quote = FALSE)
write.csv(data.g, "Data_N/data_g.csv", quote = FALSE)
write.csv(data.y, "Data_N/data_y.csv", quote = FALSE, row.names = FALSE)
write.csv(as.data.frame(POS), "Data_N/data_pos.csv", row.names = FALSE)





# data.. <- data %>% group_by(PTID) %>% filter(n() > 1) %>% mutate(diff = AGE - min(AGE)) %>% slice(2:n())
