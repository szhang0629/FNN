library(refund)
Error.flm2 <- function(Y.train, G.train, Y.test, G.test, pos, loc){
  odr. <- order(loc)
  loc <- loc[odr.]
  Y.train <- Y.train[, odr.]
  Y.test <- Y.test[, odr.]
  data <- data.frame(X1=I(G.train), Y=I(Y.train),  fix.empty.names = F)
  data.test <- data.frame(X1=I(G.test), Y=I(Y.test), fix.empty.names = F)
  # nz <- ceiling(ncol(G.train)*(1-exp(-nrow(G.train)/ncol(G.train))))
  m1 <- pffr(Y ~ ff(X1, xind = pos, splinepars =
                      list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                           # k = c(round(nz/3), round(length(loc)/2)))),
             k = c(round(20*log(ncol(G.train))/log(80)), 20))), ## !!!
             yind = loc, data = data)
  # m1 <- pffr(Y ~ ff(X1, xind = pos), yind = loc, data = data, 
  #            bs.yindex = list(bs = "ps", k = round(10), m=c(2, 1)), 
  #            bs.int = list(bs = "ps", k = 10, m=c(2, 1)))
  Y.train. <- predict(m1, data)
  Y.test. <- predict(m1, data.test)
  train <- cost.(Y.train, Y.train.)
  test <- cost.(Y.test, Y.test.)
  cor1 <- cor(c(unlist(Y.train)), c(unlist(Y.train.)))
  cor2 <- cor(c(unlist(Y.test)), c(unlist(Y.test.)))
  result <- data.frame(train = train, test = test, lambda = 0, 
                       cor1 = cor1, cor2 = cor2)
  return(result)
}