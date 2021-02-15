divide <- function(IDs, seed = 629, type = "index", ratio = 0.8, 
                   names = c("train", "test")) {
  IDs <- gsub("\\..*", "", IDs)
  ID <- unique(IDs)
  N <- length(ID)
  set.seed(seed)
  idx <- sample(1:N)
  train = list()
  test = list()
  if (type == "name") {
    train[[length(train) + 1]] <- ID[sort(idx[1:round(ratio*N)])]
    test[[length(test) + 1]] <- ID[sort(idx[(round(ratio*N) + 1):N])]
  }
  train[[length(train) + 1]] <- which(IDs %in% ID[idx[1:round(ratio*N)]])
  test[[length(test) + 1]] <- which(IDs %in% ID[idx[(round(ratio*N) + 1):N]])
  assign(names[[1]], train)##not finished
  assign(names[[2]], test)
  return(mget(names))
}
sep <- function(objects, groups){
  name1 <- names(groups)[[1]]
  name2 <- names(groups)[[2]]
  names.in.list <- list()
  name <- names(objects)
  for (i in 1:length(objects)) {
    name. <- sub("\\..*", "", name[i])
    assign(paste0(name., ".", name1), subs(objects[[i]], groups[[1]]))
    assign(paste0(name., ".", name2), subs(objects[[i]], groups[[2]]))
    names.in.list.new <- mget(c(paste0(name., ".", name1),
                                paste0(name., ".", name2)))
    names.in.list <- c(names.in.list, names.in.list.new)
  }
  return(names.in.list)
}
subs <- function(object, rowss) {
  if (is.null(nrow(object))) return(object)
  for (rows in rowss) {
    if (is.numeric(rows) || is.logical(rows)) {
        return(object[rows, , drop = FALSE])
    } else {
      if (is.data.frame(object)) {
        name <- object$PTID
        rows. <- name %in% rows
      } else {
        name <- rownames(object)
        if (!is.null(name)) {
          name. <- sub("\\..*", "", name)
          rows. <- name. %in% rows
        } else rows. <- FALSE
      }
      if (any(rows.)) return(object[rows., , drop = FALSE])
    }
  }
}
divide.cv <- function(IDs, seed = 629, n = 10) {
  IDs <- gsub("\\..*", "", IDs)
  ID <- unique(IDs)
  N <- length(ID)
  set.seed(seed)
  idx <- sample(1:N)
  groups <- list()
  for (i in 1:n) 
    groups[[i]] <- which(IDs %in% ID[idx[(round((i-1)*N/n) + 1):round(i*N/n)]])
  return(groups)
}