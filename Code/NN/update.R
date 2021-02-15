RMS <- function(g, epsilon = 1e-6){
  return(sqrt(g + epsilon))
}
adadelta <- function(output, grads.g, grads.p, rho = 0.95) {
  for (i in 1:length(grads.g)) {
    grads. <- grads.g[[i]] + grads.p[[i]]
    output$eg[[i]] <- rho * output$eg[[i]] + (1 - rho) * (grads.)^2
    dAc <- -RMS(output$ex[[i]]) / RMS(output$eg[[i]]) * grads.
    output$ex[[i]] <- rho * output$ex[[i]] + (1 - rho) * dAc^2
    output$Ac[[i]] <- output$Ac[[i]] + dAc
  }
  return(output)
}
Error.p <- function(Ac, lambda) {
  result <- rep(0, length(Ac))
  for (i in 1:length(Ac))
    result[i] <- sum(Ac[[i]]^2) * lambda
  return(result)
}
pen <- function(Ac, lambda) {
  grads <- list()
  for (i in 1:length(Ac))
    grads[[i]] <- 2 * lambda * Ac[[i]]
  return(grads)
}
# Error.p <- function(Ac, lambda) {
#   result <- rep(0, length(Ac))
#   for (i in 1:length(Ac))
#     result[i] <- sum(abs(Ac[[i]])) * lambda
#   return(result)
# }
# pen <- function(Ac, lambda) {
#   grads <- list()
#   for (i in 1:length(Ac)) 
#     grads[[i]] <- 2 * lambda * sign(Ac[[i]])
#   return(grads)
# }