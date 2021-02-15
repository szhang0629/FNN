library(tidyverse)
iplot <- function(name, levels1 = 
                    c("scalar", "vector", "matrix"),
                    # c("polyno", "logist", "linear"),
                    # c("80", "200", "500"),
                  levels2 = c("mse1", "mse2"),
                  # levels2 = c("train", "test"), 
                  levels3 = c("FLM", "NN", "FNN")) {
  # "polyno", "logist", "linear"
  
  D0 <- read.csv(paste0("Result/", name, ".csv"))
  D0[is.na(D0)] <- 0
  D1 <- D0 %>% group_by(type, noise, index, method) %>% 
    slice(which.min(train)) %>% ungroup()
  D1 <- D1[, c('type', 'noise', 'index', 'method', levels2)]
  D1 <- D1 %>% rename(train = mse1, test = mse2)
  levels2 = c("train", "test")
  D2 <- gather(D1, key = class, value = "mse", train, test)
  # D2$type <- recode(D2$type, fct_1d = "vector", fct_2d = "matrix")
  D2$noise <- recode(D2$noise, "1" = "0.3", "2" = "0.6", "3" = "1.2")
  # D2 <- gather(D1, key = class, value = "mse", train, test)
  D2$method <- recode(D2$method, FLM1 = "FLM", NN1 = "NN", FN1 = "FNN")
  D2 <- D2[D2$method %in% levels3, ]
  D2$method <- factor(D2$method, levels = levels3)
  D2$type <- factor(D2$type, levels = levels1)
  D2$class <- factor(D2$class, levels = levels2)
  D2 <- D2 %>% group_by(noise, method, type, class) %>% 
    mutate(mse = pmin(mse, quantile(mse, 0.95)))
  D2 <- D2 %>% group_by(noise, method, type, class) %>% 
    mutate(mse = pmax(mse, quantile(mse, 0.05)))
  ggplot(data = D2,  mapping = aes(x = class, y = mse, fill = method)) +
    facet_grid(noise ~ type, scales = "free_y") +
    geom_boxplot()
  # ggsave(paste0(name, ".pdf"), width = 8, height = 8)
  ggsave(file=paste0(name, ".eps"), device="eps", width = 8, height = 8)
}

iplot("Sim1", levels1 = c("scalar", "vector", "matrix"))
iplot("Sim2", levels1 = c("polyno", "logist", "linear"))
iplot("Sim3", levels1 = c("80", "200", "500"))
