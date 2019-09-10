gData <- function(index, noise = 0.1, n = 200, p = 500, loc = 0.5, fct.){
  set.seed(index)
  Data. <- Data.Input(n = n, p = p)
  # f <- Fun.Output(Data.$G, Data.$pos, index)
  fct <- Fun.Output(fct.)
  f <- fct(Data.$G, Data.$pos, index)
  Y <- Data.Output(f, loc, noise)
  Data <- list(G = Data.$G, f = f, Y = Y, pos = Data.$pos)
  return(Data)
}