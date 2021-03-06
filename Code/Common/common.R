cost. <- function(Y, Y.hat = NULL){
  return(mean((Y - Y.hat)^2))
}
prt <- function(x, path, col.names = F, pscn = F) {
  if (pscn)
    print(x)
  write.table(x, path, append = T, quote = F, sep = ",",
              row.names = F, col.names = col.names)
}
Output. <- function(x = NULL, path) {
  folders <- strsplit(path, "/")[[1]]
  ipath <- folders[[1]]
  for (i in 2:(length(folders) - 1)) {
    ipath <- paste(ipath, folders[i], sep = "/")
    if (!dir.exists(ipath))
      dir.create(ipath)
  }
  ipath <- paste(ipath, folders[length(folders)], sep = "/")
  if (!file.exists(ipath)) {
    if (is.null(x)) file.create(ipath)
    else prt(x, ipath)
  }
}
rep.row <- function(x,n){
  matrix(rep(x,each = n),nrow = n)
}
rep.col <- function(x,n){
  matrix(rep(x,each = n), ncol = n, byrow = TRUE)
}
Error <- function(Y.train, Y.test, Y.train., Y.test.) {
  if (is.matrix(Y.train)) {
    val <- colMeans(Y.train)
    val1 <- cost.(Y.train, rep.row(val, nrow(Y.train)))
    val2 <- cost.(Y.test, rep.row(val, nrow(Y.test)))
  } else {
    val <- mean(Y.train)
    val1 <- cost.(Y.train, val)
    val2 <- cost.(Y.test, val)
  }
  mse1 <- cost.(c(Y.train), c(Y.train.))
  mse2 <- cost.(c(Y.test), c(Y.test.))
  train <- mse1/val1
  test <- mse2/val2
  return(data.frame(train = train, test = test, mse1 = mse1, mse2 = mse2))
}
cbb <- function(norder, pos, nbasis = NA) {
  pos <- sort(pos)
  breaks. <- ceiling(seq(1, length(pos), length.out = nbasis - norder + 1))
  breaks. <- pos[breaks.]
  breaks <- (c(0, breaks.) + c(breaks., 1))/2
  return(create.bspline.basis(rangeval = 0:1, breaks = breaks, norder = norder))
}
seq. <- function(from, to, len) {
  return(10^(seq(from, to, length.out = len + 2)))
}
func.plot <- function(coef, basisobj) {
  x <- (0:99)/100 + 0.005
  plot(x, eval.basis(x, basisobj) %*% coef)
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