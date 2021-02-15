library(tidyverse)
iplot <- function(name) {
  .vari <- c('idx', 'method', 'train', 'test')
  # .vari <- c('idx', 'method', 'mse1', 'mse2')
  data <- read.csv(paste0("Result/", name, ".csv"))[, .vari]
  data <- data[data$method != "Base",]
  data <- data %>% group_by(idx) %>% filter(n() == 5) %>% ungroup()
  data <- data %>% gather("class", "mse", train, test)
  # data <- data %>% gather("class", "mse", mse1, mse2)
  data$class <- recode(data$class, mse1 = "train", mse2 = "test")
  data$class <- factor(data$class, levels = c("train", "test"))
  # data <- data %>% group_by(method) %>% 
  #   mutate(train = pmin(train, quantile(mse, 0.975))) %>% 
  #   mutate(train = pmax(train, quantile(mse, 0.025))) %>% 
  #   mutate(test = pmin(test, quantile(mse, 0.975))) %>% 
  #   mutate(test = pmax(test, quantile(mse, 0.025)))
  data$method <- recode(data$method, NN1 = "NN1", FN1 = "FNN1", 
                      NN2 = "NN2", FN2 = "FNN2", FLM1 = "FLM")
  data$method <- factor(data$method, levels = 
                          c("FLM", "NN1", "NN2", "FNN1", "FNN2"))
  levels(data$class) <- c("train", "test")
  ggplot(data = data,  mapping = aes(x = class, y = mse, fill = method)) +
    # coord_cartesian(ylim = c(0.85, 1.1)) +
    # labs(title = paste0("MSE of ", substr(name, 1, nchar(name) - 1), 
    #                     " on Smoking Frequency")) +
    geom_boxplot()

  ggsave(paste0(name, ".eps"), width = 6, height = 4, device="eps")
}

iplot("Real1")
iplot("Real2")
