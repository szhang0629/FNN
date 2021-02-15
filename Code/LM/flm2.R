library(refund)
Error.flm2 <- function(Y.train, G.train, Y.test, G.test, pos, loc){
  odr. <- order(loc)
  loc <- loc[odr.]
  Y.train <- Y.train[, odr.]
  Y.test <- Y.test[, odr.]
  data <- data.frame(X1 = I(G.train), Y = I(Y.train),  fix.empty.names = F)
  data.test <- data.frame(X1 = I(G.test), Y = I(Y.test), fix.empty.names = F)
  # nz <- ceiling(ncol(G.train)*(1-exp(-nrow(G.train)/ncol(G.train))))
  m1 <- pffr(Y ~ ff(X1, xind = pos
                    ,splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                    k = c(20, 20))
             ), yind = loc, data = data)
  # m1 <- pffr(Y ~ ff(X1, xind = pos), yind = loc, data = data, 
  #            bs.yindex = list(bs = "ps", k = round(10), m=c(2, 1)), 
  #            bs.int = list(bs = "ps", k = 10, m=c(2, 1)))
  Y.train. <- predict(m1, data)
  Y.test. <- predict(m1, data.test)
  error <- Error(Y.train, Y.test, Y.train., Y.test.)
  return(error)
}