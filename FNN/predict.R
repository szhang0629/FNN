pred. <- function(Ac, X, G, A, Bases, int, loc){ 
  layers <- forw.prop(Ac, X, G, Bases, A, int, loc)
  return(layers[[length(layers)]]$G)
}
pred <- function(Ac, X, G, A, Bases = NULL, pos = NULL, loc = NULL, int = NULL){
  if (!is.null(pos)) {
    int <- int.fun(Bases, int)
    idy <- idx.names(subset(loc, select = -PTID))
    Bases <- Bf2m(Bases, pos, subs(subset(loc, select = -PTID), 
                                   match(1:max(idy), idy)))
    idx <- idx.names(loc$PTID, rownames(G))
    loc <- data.frame(idx = idx, idy = idy)
  }
  return(pred.(Ac, X, G, A, Bases, int = int, loc))
}