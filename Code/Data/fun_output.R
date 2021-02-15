Fun.Output <- function(fct.name = "polyno") {
  fct <- switch(fct.name, polyno = .polyno, logist = .logist, linear = .linear)
  force(fct)
  return(fct)
}