pred. <- function(Ac, X, G, A, Bases = NULL, loc){ 
  layers <- forw.prop(Ac, X, G, Bases, A, loc)
  return(layers[[length(layers)]]$G)
}
pred <- function(Ac, X, G, A, Bases = NULL, pos = NULL, loc = NULL){
  if (is.null(Bases)) {
    Bases <- Bf2m(rep(list(.bb), length(A) + 1), pos, 
                  subset(loc, select = -(PTID)))
    loc$PTID <- idx.names(loc$PTID, rownames(G))
    loc$loc <- idx.names(subset(loc, select = -PTID))
    loc <- subset(loc, select = c(PTID, loc))
  }
  # # Bases <- Bf2m(Bases, pos, subset(loc, select = -(PTID)))
  # if (!is.numeric(loc$PTID)) {
  #   loc$PTID <- idx.names(loc$PTID, rownames(G))
  #   loc$loc <- idx.names(subset(loc, select = -(PTID)))
  #   loc <- subset(loc, select = c(PTID, loc))
  # }
  return(pred.(Ac, X, G, A, Bases, loc))
}