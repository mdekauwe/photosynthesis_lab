#' General utility functions for the photosynthesis model
#'
#' These functions support core model calculations, including math utilities
#' like approximate equality checks and root solving.
#'
#' @author Martin De Kauwe
#' @email mdekauwe@gmail.com
#' @date 20/06/2025
#' @name utils
#' @keywords internal
#'
NULL


#' Check if two numeric values are approximately equal
#'
#' @param a Numeric value
#' @param b Numeric value
#' @param rel_tol Relative tolerance (default 1e-9)
#' @param abs_tol Absolute tolerance (default 0.0)
#'
#' @return Logical. TRUE if values are considered close, FALSE otherwise.
is_close <- function(a, b, rel_tol=1e-09, abs_tol=0.0) {
  return ( abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol) )
}

#' Solve a quadratic equation and return one real root
#'
#' @param a Coefficient of x^2
#' @param b Coefficient of x
#' @param c Constant term
#' @param large Logical. If TRUE, returns the larger root; otherwise, the 
#'              smaller root.
#'              
#' @return Numeric. One real root of the quadratic equation.
quadratic <- function(a, b, c, large=FALSE) {
 
  # discriminant
  d <- b**2.0 - 4.0 * a * c

  if (d < 0) stop("Imaginary root found")
  
  if (is_close(a, 0)) {
    if (!is_close(b, 0)) {
      return (-c / b)
    } else if (is_close(c, 0)) {
      return (0.0)  # All coefficients zero, infinite roots 
    } else {
      stop("Cannot solve: no valid root")
    }
  }
  
  sqrt_d <- sqrt(d)
  if (isTRUE(large)) {
    root <- (-b + sqrt_d) / (2 * a)
  } else {
    root <- (-b - sqrt_d) / (2 * a)
  }
  
  return(root)
}
