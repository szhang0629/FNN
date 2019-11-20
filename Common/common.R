cost. <- function(Y, Y.hat = NULL){
  return(mean((Y - Y.hat)^2))
}
prt <- function(x, path, col.names = F) {
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
  if (is.data.frame(Y.train)) {
    Y.train <- Y.train$Y
    Y.test <- Y.test$Y
  }
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
  cor1 <- cor(c(Y.train), c(Y.train.))
  cor2 <- cor(c(Y.test), c(Y.test.))
  return(data.frame(train = train, test = test, mse1 = mse1, mse2 = mse2, 
                    cor1 = cor1, cor2 = cor2))
}