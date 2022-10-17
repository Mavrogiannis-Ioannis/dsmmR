#' lambda genome
#'
#' Complete data from phage genome [WT71] of length 48502. (needs a citation...)
#'
#' @name lambda
#' @docType data
#'
#' @usage data(lambda)
#' data("lambda")
#' # With the following we do not load the `dsmmR` package.
#' data(lambda, package = "dsmmR")
#' data("lambda", package = "dsmmR")
#'
#' @format A vector object of type \code{"character"} and length of 48502.
#' It has class of \code{"Rdata"}.
#'
#' @keywords datasets
#'
#' @seealso \code{\link[utils]{data}}
#'
#' @examples
#' data(lambda, package = "dsmmR")
#' class(lambda)
#' # > "SeqFastadna"
#' sequence <- c(lambda) # Convert to "character" class
#' str(sequence)
#' # > chr [1:48502] "g" "g" "g" "c" "g" "g" "c" "g" "a" "c" "c" "t" ...
NULL

