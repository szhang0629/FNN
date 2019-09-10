FNN. <- function(Y, E, output, Bases, A, A.prime, lambda = 0, P = NULL, 
                 int = T) {
  if (is.list(Y)) lambda <- lambda/length(Y)
  else lambda <- lambda <- lambda/nrow(Y)
  layers <- forw.prop(E, output$Ac, Bases, A, int = int)
  grads <- back.prop(Y, output$Ac, layers, Bases, A, A.prime, int = int)
  output <- adadelta(output, grads, lambda, P)
  return(output)
}
FNN <- function(Y, E, loc = NULL, pos = NULL, Bases., A, lambda. = 0, 
                Ac.init = NULL){
  P <- pen.Bases(Bases.)
  int <- int.fun(Bases.)
  Bases <- Bf2m(Bases., pos, loc)
  Ac <- init(Bases, Ac.init, seed = 1, X.col = ncol(E$X))
  e0 <- init(Bases, X.col = ncol(E$X), para.ub = 0)
  output <- list(Ac = Ac, eg. = e0, ex. = e0)
  A.prime <- list()
  for (i in 1:length(A))
    A.prime[[i]] <- Deriv(A[[i]])
  return(FNN.p(Y, E, output, Bases, A, A.prime, lambda., P, int))
}