
##' Cubic spline interpolation between a set of points with given
##' derivatives at those points
##'
##' Extrapolation using a constant function from the boundary value.
##' This keeps the gradient continuously 0 at and beyond the boundary
##'
##' @details From stats:::splinefunH0 TODO full credits
##'
##' @param x Vector of points at which to evaluate the interpolated function
##' @param x0 Vector of x-positions of points to be interpolated
##' @param y0 Vector of y-positions of points to be interpolate
##' @param m Vector of derivatives at the points (`x0`, `y0`)
##' @param lower Values around here will be rounded
##' @param upper Values around here will be rounded
##' @noRd
hermite <- function(x, x0, y0, m, lower=0, upper=1){
  dx <- x0[-1L] - x0[-length(x0)]
  i <- findInterval_soft(x, x0)
  h <- dx[i]
  t <- (x - x0[i])/h
  t1 <- t - 1
  h01 <- t*t*(3 - 2*t)
  h00 <- 1 - h01
  tt1 <- t*t1
  h10 <- tt1*t1
  h11 <- tt1*t
  ret <- y0[i]*h00 + h*m[i]*h10 + y0[i+1]*h01 + h*m[i+1]*h11
  ret[ret<lower & lower-ret < .Machine$double.eps] <- lower
  ret[ret>upper & ret-upper < .Machine$double.eps] <- upper
  ret[x < x0[1]] <- y0[1]
  ret[x > x0[length(x0)]] <- y0[length(x0)]
  ret
}

findInterval_soft <- function(x, vec){
  i <- findInterval(x, vec, rightmost.closed=TRUE, all.inside=TRUE)
#  i[((x < vec[1]) &  (x > vec[1] - .Machine$double.eps))] <- 1
#  i[((x > vec[length(vec)]) & (x < vec[length(vec)] + .Machine$double.eps))] <- length(vec)-1
  i
}

##' Pointwise derivatives required to fit a Hermite spline to interpolate data
##' Based on observed second order differences in the data (or first order at the ends)
##' @noRd
hermite_point_derivs <- function(x, y, lower=0, upper=1, zero_threshold=1e-01){
  ## TODO arg checks
  n <- length(y)
  m_lower <- (y[2] - y[1]) / (x[2] - x[1])
  m_upper <- (y[n] - y[n-1]) / (x[n] - x[n-1])
  m_mid <- diff(y,lag=2)/diff(x,lag=2)
  m <- c(m_lower, m_mid, m_upper)
  m[all_equal_vec(y, lower)] <- 0
  m[all_equal_vec(y, upper)] <- 0
  m[y<zero_threshold] <- 0
  m
}
##' @param target vector
##' @param current scalar
##' @noRd
all_equal_vec <- function(target, current){
  Vectorize(function(x, y)isTRUE(all.equal(x,y)))(target,current)
}

##' Derivative of spline interpolation function
##' @inheritParams hermite
##' @noRd
Dhermite <- function(x, x0, y0, m, lower=0, upper=1){
  dx <- x0[-1L] - x0[-length(x0)]
  i <- findInterval_soft(x, x0)
  h <- dx[i]
  t <- (x - x0[i])/h
  t1 <- t - 1
  h01 <- -6 * t * t1
  h10 <- (3 * t - 1) * t1
  h11 <- (3 * t - 2) * t
  ret <- (y0[i+1] - y0[i])/h*h01 + m[i]*h10 + m[i+1]*h11
  ret[x < x0[1]] <- 0
  ret[x > x0[length(x0)]] <- 0
  ret
}
