diffrac <- function (x, y = NULL) 
{
  
  if (is.data.frame(y)) 
    y <- as.matrix(y)
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.matrix(x) && is.null(y)) 
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y))) 
      stop("'y' must be numeric")
    stopifnot(is.atomic(y))
  }
 
  Rank <- function(u) {
    ## take care not to drop dims on a 0- or 1-row matrix
    if(length(u) == 0L) u
    else if(is.matrix(u)) {
      if(nrow(u) > 1L) apply(u, 2L, rank, na.last="keep") else row(u)
    } else rank(u, na.last="keep")
  }

  if (is.null(y)){  # rank correlations and "pairwise.complete.obs"; the hard case
    ## Based on contribution from Shigenobu Aoki.
    ## matrix
      ncy <- ncx <- ncol(x)
      if (ncx == 0) 
        stop("'x' is empty")
      r <- matrix(0, nrow = ncx, ncol = ncy)
      for (i in seq_len(ncx)) {
        for (j in seq_len(i)) {

          x2 <- x[, i]
          y2 <- x[, j]
          ok <- complete.cases(x2, y2)

          r[i, j] <-  if(any(ok)) sum(abs(x2 - y2), na.rm = TRUE) else NA


        }


      }
      
      r <- r + t(r) - diag(diag(r)) # Make symmetric 
      rownames(r) <- colnames(x)
      colnames(r) <- colnames(x)
      
      return(r)

  }
    ## vector/matrix x vector/matrix
    else {
      if (length(x) == 0L || length(y) == 0L) 
        stop("both 'x' and 'y' must be non-empty")
      matrix_result <- is.matrix(x) || is.matrix(y)
      if (!is.matrix(x)) 
        x <- matrix(x, ncol = 1L)
      if (!is.matrix(y)) 
        y <- matrix(y, ncol = 1L)
      ncx <- ncol(x)
      ncy <- ncol(y)
      r <- matrix(0, nrow = ncx, ncol = ncy)
      for (i in seq_len(ncx)) {
        for (j in seq_len(ncy)) {
          x2 <- x[, i]
          y2 <- y[, j]
          ok <- complete.cases(x2, y2)
          r[i, j] <-  if(any(ok)) sum(abs(x2 - y2), na.rm = TRUE) else NA

      }
      rownames(r) <- colnames(x)
      colnames(r) <- colnames(y)
      if (matrix_result) 
        r
      else drop(r)
      }
  return(r)
      }
}
