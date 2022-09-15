#' Logit link function
#'
#' @param x value of x to be transformed, x e (0, 1)
#'
#' @return logit value of x, x unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(x, logit(x), type="l")
logit <- function(x) {
  return(qlogis(x))
}

#' Inverse-Logit link function, equal to Logistic link
#'
#' @param x value of x to be transformed, x unbound
#'
#' @return inverse logit value of x, result is e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -5, to = 5, length.out = 100)
#' plot(x, inv_logit(x), type="l")
inv_logit <- function(x) {
  return(plogis(x))
}

#' Logistic link function, equivalent to Inverse-Logit link
#'
#' @param x value of x to be transformed, unbound
#'
#' @return logistic value of x, result ise (0, 1)
#' @export
#'
#' @examples x <- seq(from = -100, to = 100, length.out = 100)
#' plot(x, logistic(x), type="l")
logistic <- function(x) {
  return(inv_logit(x))
}


#' Complementary-Log-Log link function
#'
#' @param x value of x to be transformed, x e (0, 1)
#'
#' @return cloglog value of x, result is unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(x, cloglog(x), type="l")
cloglog <- function(x) {
  log(-log1p(-x))
}


#' Inverse CLogLog link function
#'
#' @param x value of x to be transformed, x unbound
#'
#' @return inverse-cloglog value of x, result is e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -3, to = 1, length.out = 100)
#' plot(x, inv_cloglog(x), type="l")
inv_cloglog <- function(x) {
  return(1 - exp(-exp(x)))
}


#' Cauchit link function
#'
#' @param x value of x to be transformed, not defined for x of whole integer
#'
#' @return cauchit value of x, result is unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(x, cauchit(x), type="l")
cauchit <- function(x) {
  return(qcauchy(x))
}


#' Inverse Cauchit link function equivalent to pcauchy
#'
#' @param x value of x to be transformed, any real scalar or vector allowed
#'
#' @return inverse cauchit of x, result is e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -10, to = 10, length.out = 100)
#' plot(x, inv_cauchit(x), type="l")
inv_cauchit <- function(x) {
  return(pcauchy(x))
}

#' Gaussion Error function
#'
#' @param x value to be transformed, x unbound
#'
#' @return erf function of x, result e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -2, to = 2, length.out = 100)
#' plot(x, erf(x), type="l")
erf <- function(x) {
  return(2 * pnorm(x * sqrt(2)) - 1)
}


#' Softplus link function
#'
#' @param x value to be transformed, x positive unbound
#'
#' @return softplus function of x, result is unbound
#' @export
#'
#' @examples x <- seq(from = -5, to = 5, length.out = 100)
#' plot(x, softplus(x), type="l")
softplus <- function(x) {
  return(log(exp(x) - 1))
}


#' Inverse Softplus link function
#'
#' @param x value to be transformed, x is unbound
#'
#' @return inv_softplus of x, result is positive unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 5, length.out = 100)
#' plot(x, softplus(x), type="l")
inv_softplus <- function(x) {
  return(log(exp(x) + 1))
}


#' Numeric vector check
#'
#' @param num Numeric vector to be checked
#' @param len Length of vector, default argument is 1
#'
#' @return Boolean, whether num was numeric and of correct size
#' @export
#'
#' @examples isNum_len(c(1.1, 2.2), 2) # should be TRUE
#' isNum_len(0.2) # should be TRUE
#' isNum_len(0.2, 2) # should be FALSE, wrong length
#' isNum_len("r") # should be FALSE, not numeric
#' isNum_len(c("r", 0.2), 2) # should be FALSE, partially not numeric
#' isNum_len(c("r", 0.2), 3) # should be FALSE, both not numeric and wrong length
#'
isNum_len <- function(num, len=1) {
  value <- (!any(is.na(num))) && is.numeric(num) && length(num) == len
  return(isTRUE(value))
}

#' Integer vector check
#'
#' @param int Integer vector to be checked
#' @param len Length of vector, default argument is 1
#'
#' @return Boolean, whether int was Integer and of correct size
#' @export
#'
#' @examples isInt_len(c(1, 2), 2) # should be TRUE
#' isInt_len(1) # should be TRUE
#' isInt_len(1, 2) # should be FALSE, wrong length
#' isInt_len("r") # should be FALSE, not integer
#' isInt_len(0.2) # should be FALSE, not integer
#' isInt_len(c("r", 1), 2) # should be FALSE, partially not integer
#' isInt_len(c("r", 1), 3) # should be FALSE, both not integer and wrong length
isInt_len <- function(int, len=1) {
  if(isTRUE((!any(is.na(int))) && is.numeric(int))) {
    # only check for integer, if the type is numeric!
    value <- all(int %% 1 == 0) && length(int) == len
    return(isTRUE(value))
  }
  else {
    return(FALSE)
  }
  # One might also do this all in a single AND beginning with is.numeric.
  # In a single test, this worked fine, given if !is.numeric, the other boolean,
  # checks were not done (because FALSE & x <=> FALSE)
  # this did prevent errors (from "r" %% 1 == 0) at the least

  # But I would prefer only checking numerics,
  # given logical AND may not necassarily follow this behaviour!
}

#' Check, if the input is a single string
#'
#' @param input String argument
#'
#' @return Is a string and only one string
#' @export
#'
#' @examples isSingleString("abc")  # should be TRUE
#' isSingleString(c("abc")) # should be TRUE
#' isSingleString(c("abc", "def")) # should be FALSE, not a single string
#' isSingleString(1) # should be FALSE, not a string
#' isSingleString(c("abc", 1)) # should be FALSE, partially not a string and
#' # also wrong length
isSingleString <- function(input) {
  value <- (!any(is.na(input))) && is.character(input) && length(input) == 1
  return (isTRUE(value))
}

#' Data limit function
#'
#' @param data Data to be limited
#' @param limits Limits to be used. Vector with 2 real entries, limits[1] <= limits[2]
#' If the lower bound does not have to be restricted, set it to NA and vice versa.
#' Sets data outside those bounds to those bounds.
#'
#' @return data limited by the limits
#' @export
#'
#' @examples input <- c(1, 2, 3, 4)
#' print(input)
#' print(limit_data(input, c(2, 3)))
#' print(limit_data(input, c(2, NA)))
limit_data <- function(data, limits) {
  # check that the limit is usable
  if(length(limits) != 2) {
    stop("If the limits is to be used, it has to be of size 2.")
  }
  if(!isNum_len(data, len=length(data))) {
    stop("Some data was not numeric, or was NA")
  }

  # if so, use the applicable limit (If one uses to na, well. What are you trying to achieve? :)
  if(isNum_len(limits, 2)) {
    if(limits[1] > limits[2]) {
      stop("In limit_data, the first limit is the lower limit, so it has to be
           smaller than the second limit.")
    }
  }

  # isNum_len will certailny return false, if NA
  if(isNum_len(limits[1])) {
    data[data < limits[1]] <- limits[1]
  }
  if(isNum_len(limits[2])) {
    data[data > limits[2]] <- limits[2]
  }

  # now return the data
  return(data)
}
