#' lambda genome
#'
#' Staphylococcus aureus is a widespread and versatile bacterium
#' that can infect humans as well as animals.
#' It can also lead to some life-threatening diseases like pneumonia.
#' The following is the complete DNA sequence
#' from an isolate of the S. aureus bacterium, with number [WT71],
#' that was considered on a study on diseased Eurasian Beavers.
#'
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
#' @seealso \code{\link[utils:data]{data}}
#' @references S. Monecke & A. Feßler & S.Eurasian Beavers ( Burgold-Voigt &
#'  H. Krüger-Haker & K. Muehldorfer & G. Wibbelt &
#'  E. Liebler-Tenorio & M. Reinicke & S. Braun &
#'  D. Hanke & C. Diezel & E. Müller & I. Loncaric &
#'  S. Schwarz & R. Ehricht (2021).
#'    Staphylococcus aureus isolates from Eurasian Beavers (Castor fiber) carry a novel
#'    phage-borne bicomponent leukocidin related to the Panton-Valentine leukocidin.
#'    Scientific Reports. 11. 10.1038/s41598-021-03823-6.
#' @examples
#' data(lambda, package = "dsmmR")
#' class(lambda)
#' sequence <- c(lambda) # Convert to "character" class
#' str(sequence)
NULL

