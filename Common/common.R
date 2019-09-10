cost. <- function(Y, Y.hat = NULL){
  if (is.list(Y)) {
    Y <- unlist(Y)
    Y.hat <- unlist(Y.hat)
  }
  if (is.null(Y.hat))
    return(mean(Y^2))
  else
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