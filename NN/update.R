RMS <- function(g, epsilon = 1e-6){
  return(sqrt(g + epsilon))
}
adadelta <- function(output, grads, lambda = 0, P = NULL, P0 = NULL, 
                     rho = 0.95) {
  eg <- list()
  ex <- list()
  Ac <- output$Ac
  eg. <- output$eg.
  ex. <- output$ex.
  for (i in 1:length(grads)) {
    grads. <- grads[[i]] + lambda * Ac[[i]]
    eg[[i]] <- rho * eg.[[i]] + (1 - rho) * (grads.)^2
    dAc <- - RMS(ex.[[i]]) / RMS(eg[[i]]) * grads.
    ex[[i]] <- rho * ex.[[i]] + (1 - rho) * dAc^2
    Ac[[i]] <- Ac[[i]] + dAc
  }
  return(list(Ac = Ac, eg. = eg, ex. = ex))
}
integ <- function(D, B, int = F){
  if (is.null(D))
    return(NULL)
  if (int) return(as.matrix(D) %*% B / nrow(B))
  else return(as.matrix(D) %*% B)
}
sigmoid <- function(x){
  (1 + exp(-x))^(-1)
}
linear <- function(x){
  x
}
relu <- function(x){
  (abs(x) + x) / 2
}