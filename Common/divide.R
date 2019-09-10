divide. <- function(N, seed = 629, ratio = 0.8, names = c("train", "test")){
  set.seed(seed)
  idx <- sample(1:N)
  assign(names[[1]], sort(idx[1:ceiling(ratio*N)]))
  assign(names[[2]], sort(idx[(ceiling(ratio*N) + 1):N]))
  return(mget(names))
}
divide.g <- function(IDs, seed = 629, ratio = 0.8, names = c("train", "test")) {
  ID <- unique(IDs)
  divide <- divide.(length(ID), seed)
  .train <- divide[[1]]
  .test <- divide[[2]]
  train <- grep(paste(ID[.train], collapse = "|"), IDs)
  test <- grep(paste(ID[.test], collapse = "|"), IDs)
  assign(names[[1]], sort(train))
  assign(names[[2]], sort(test))
  return(mget(names))
}
divide <- function(x, seed = 629, ratio = 0.8, names = c("train", "test")) {
  if (is.list(x)) return(divide.(length(x), seed, ratio, names))
  else if (min(nchar(rownames(x))) > 10 && !is.null(rownames(x)))
    return(divide.g(substr(rownames(x), 1, 10), seed, ratio, names))
  else return(divide.(nrow(x), seed, ratio, names))
}
split <- function(objects, groups){
  n <- length(groups[[1]]) + length(groups[[2]])
  name1 <- names(groups)[[1]]
  name2 <- names(groups)[[2]]
  names.in.list <- list()
  name <- names(objects)
  for (i in 1:length(objects)) {
    name. <- sub("\\..*", "", name[i])
    assign(paste0(name., ".", name1), subs(objects[[i]], groups[[1]], n))
    assign(paste0(name., ".", name2), subs(objects[[i]], groups[[2]], n))
    names.in.list.new <- mget(c(paste0(name., ".", name1),
                                paste0(name., ".", name2)))
    names.in.list <- c(names.in.list, names.in.list.new)
  }
  return(names.in.list)
}
subs <- function(object, rows, n, tolerance = 1) {
  if (is.matrix(object) || is.data.frame(object))
    if (nrow(object) == n)
      return(object[rows, , drop = FALSE])
  if (!is.basis(object) && is.list(object)) {
    if (length(object) == n)
      return(object[rows])
    else if (tolerance > 0) {
      for (i in 1:length(object))
        if (!is.null(object[[i]]))
          object[[i]] <- subs(object[[i]], rows, n, tolerance - 1)
        return(object)
    }
  }
  return(object)
}