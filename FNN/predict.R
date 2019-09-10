pred. <- function(Ac, E, A, Bases = NULL, int = T){ 
  layers <- forw.prop(E, Ac, Bases, A, int = int)
  result <- layers[[length(layers)]]$E
  return(result)
}
pred <- function(Ac, E, A, Bases = NULL, pos = NULL, loc = NULL, int = NULL){
  int <- int.fun(Bases, int)
  Bases <- Bf2m(Bases, pos, loc)
  if (is.list(Bases[[length(Bases)]])) {
    result <- list()
    Bases.z <- Bases[[length(Bases)]]
    N <- length(Bases.z)
    int <- int.fun(Bases, int)
    for (i in 1:N) {
      Bases[[length(Bases)]] <- Bases.z[[i]]
      result[[i]] <- pred.(Ac, subs(E, i, N), A, Bases, int)
    }
  } else 
    result <- pred.(Ac, E, A, Bases, int = int)
  return(result)
}