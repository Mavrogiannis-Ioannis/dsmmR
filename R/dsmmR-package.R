#' @title dsmmR : Estimation and Simulation of Drifting Semi-Markov Models
#'
#' @description Performs parametric and non-parametric estimation and simulation
#' of drifting semi-Markov processes. The definition of parametric and
#' non-parametric model specifications is also possible. Furthermore, three
#' different types of drifting semi-Markov models are considered. These models
#' differ in the number of transition matrices and sojourn time distributions
#' used for the computation of a number of semi-Markov kernels, which in turn
#' characterize the drifting semi-Markov kernel.
#'
#' @details
#'
#' \strong{Introduction}
#'
#' The difference between the Markov models and the semi-Markov models
#' concerns the modelling of the sojourn time distributions.
#' The Markov models (in discrete time) are modelled by a sojourn time
#' following the Geometric distribution. The semi-Markov models
#' are able to have a sojourn time distribution of arbitrary shape.
#' The further difference with a \emph{drifting} semi-Markov model,
#' is that we have \eqn{d + 1} (arbitrary) sojourn time distributions
#' and \eqn{d + 1} transition matrices (Model 1),
#' where \eqn{d} is defined as the polynomial degree.
#' Through them, we compute \eqn{d + 1} semi-Markov kernels.
#' In this work, we also consider the possibility for obtaining these
#' semi-Markov kernels with \eqn{d + 1} transition matrices and \eqn{1}
#' sojourn time distribution (Model 2) or \eqn{d + 1} sojourn time
#' distributions and \eqn{1} transition matrix (Model 3).
#'
#' \strong{Definition}
#'
#' Drifting semi-Markov processes are particular non-homogeneous semi-Markov
#' chains for which the drifting semi-Markov kernel
#' \eqn{q_{\frac{t}{n}}(u,v,l)} is defined as
#' the probability that, given at the instance \eqn{t}
#' the previous state is \eqn{u}, the next state state \eqn{v} will be
#' reached with a sojourn time of \eqn{l}:
#' \deqn{q_{\frac{t}{n}}(u,v,l) = P(J_{t}=v,X_{t}=l|J_{t-1}=u),}
#' where \eqn{n} is the model size, defined as the length of the
#' embedded Markov chain \eqn{(J_{t})_{t\in \{0,\dots,n\}}} minus the
#' last state, where \eqn{J_{t}} is the state at the instant \eqn{t} and
#' \eqn{X_{t}=S_{t}-S_{t-1}} is the sojourn time of the state \eqn{J_{t-1}}.
#'
#' The drifting semi-Markov kernel \eqn{q_{\frac{t}{n}}}
#' is a linear combination of the product of \eqn{d + 1} semi-Markov kernels
#' \eqn{q_{\frac{i}{d}}}, where every semi-Markov kernel is the product of
#' a transition matrix \eqn{p} and a sojourn time distribution
#' \eqn{f}. We define the situation when both \eqn{p} and
#' \eqn{f} are "drifting" between \eqn{d + 1} fixed points of the model
#' as Model 1, and thus we will use the exponential \eqn{(1)} as a way to
#' refer to the drifting semi-Markov kernel
#' \eqn{q_{\frac{t}{n}}^{\ (1)}} and corresponding
#' semi-Markov kernels \eqn{q_{\frac{i}{d}}^{\ (1)}} in this case.
#' For Model 2, we allow the transition matrix \eqn{p} to drift
#' but not the sojourn time distributions \eqn{f}, and for Model 3 we allow
#' the sojourn time distributions \eqn{f} to drift but not the transition
#' matrix \eqn{p}.
#' The exponential \eqn{(2)} or \eqn{(3)} will be used for signifying
#' Model 2 or Model 3, respectively.
#' In the general case an exponential will not be used.
#'
#'
#' \strong{\emph{Model 1}}
#'
#'
#' Both \eqn{p} and \eqn{f} are drifting in this case.
#' Thus, the drifting semi-Markov kernel \eqn{q_{\frac{t}{n}}^{\ (1)}} is a
#' linear combination of the product of \eqn{d + 1} semi-Markov kernels
#' \eqn{q_{\frac{i}{d}}^{\ (1)}}, which are given by:
#' \deqn{q_{\frac{i}{d}}^{\ (1)}(u,v,l)=
#' {p_{\frac{i}{d}}(u,v)}{f_{\frac{i}{d}}(u,v,l)},}
#' where for \eqn{i = 0,\dots,d} we have \eqn{d + 1} Markov transition matrices
#' \eqn{p_{\frac{i}{d}}(u,v)}
#' of the embedded Markov chain \eqn{(J_{t})_{t\in \{0,\dots,n\}}},
#' and \eqn{d + 1} sojourn time distributions
#' \eqn{f_{\frac{i}{d}}(u,v,l)}. Therefore, the drifting semi-Markov kernel
#' is described as:
#' \deqn{q_{\frac{t}{n}}^{\ (1)}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}(t)\ q_{\frac{i}{d}}^{\ (1)}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}(t)\ p_{\frac{i}{d}}(u,v)f_{\frac{i}{d}}(u,v,l),}
#' where \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree
#' \eqn{d}, which satisfy the conditions:
#' \deqn{\sum_{i=0}^{d}A_{i}(t) = 1,}
#' \deqn{A_i \left(\frac{nj}{d} \right)= 1_{\{i=j\}},}
#' where the indicator function \eqn{1_{\{i=j\}} = 1},
#' if \eqn{i = j}, \eqn{0} otherwise.
#'
#' \strong{\emph{Model 2}}
#'
#' In this case, \eqn{p} is drifting and \eqn{f} \strong{is not drifting}.
#' Therefore, the drifting semi-Markov kernel is now described as:
#' \deqn{q_{\frac{t}{n}}^{\ (2)}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}(t)\ q_{\frac{i}{d}}^{\ (2)}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}(t)\ p_{\frac{i}{d}}(u,v)f(u,v,l).}
#'
#' \strong{\emph{Model 3}}
#'
#' In this case, \eqn{f} is drifting and \eqn{p} \strong{is not drifting}.
#' Therefore, the drifting semi-Markov Kernel is now described as:
#' \deqn{q_{\frac{t}{n}}^{\ (3)}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}(t)\ q_{\frac{i}{d}}^{\ (3)}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}(t)\ p(u,v)f_{\frac{i}{d}}(u,v,l).}
#'
#'
#' \strong{Parametric and non-parametric model specifications}
#'
#' In this package, we can define parametric and non-parametric drifting
#' semi-Markov models.
#'
#' For the \emph{parametric} case, several discrete distributions are
#' considered for the modelling of the sojourn times:
#' Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial.
#' This is done from the function
#' \code{parametric_dsmm} which returns an object of the
#' S3 class (\code{dsmm_parametric}, \code{dsmm}).
#'
#' The \emph{non-parametric} model specification concerns the sojourn
#' time distributions when no assumptions are done about the
#' shape of the distributions. This is done through the function called
#' \code{nonparametric_dsmm()}, that returns an object of class
#' (\code{dsmm_nonparametric}, \code{dsmm}).
#'
#' It is also possible to proceed with a parametric or non-parametric
#' estimation for a model on an existing sequence through the function
#' \code{fit_dsmm()}, which returns an object with the S3 class
#' (\code{dsmm_fit_parametric}, \code{dsmm}) or
#' (\code{dsmm_fit_nonparametric}, \code{dsmm}) respectively, depending
#' on the given argument \code{estimation = "parametric"} or
#' \code{estimation = "nonparametric"} .
#'
#' Therefore, the \code{dsmm} class acts like a wrapper class
#' for drifting semi-Markov model specifications, while the classes
#' \code{dsmm_fit_parametric}, \code{dsmm_fit_nonparametric},
#' \code{dsmm_parametric} and \code{dsmm_nonparametric}
#' are exclusive to the functions that create the corresponding models,
#' and inherit methods from the \code{dsmm} class.
#'
#' In summary, based on an \code{dsmm} object
#' it is possible to use the following methods:
#' \itemize{
#'   \item Simulate a sequence through the function \code{simulate.dsmm()}.
#'   \item Get the drifting semi-Markov kernel
#'    \eqn{q_{\frac{t}{n}}(u,v,l)}, for any choice of \eqn{u,v,l} or \eqn{t},
#'    through the function \code{get_kernel()}.
#' }
#'
#'
#' \strong{Restrictions}
#'
#' The following restrictions must be satisfied for every drifting semi-Markov model:
#' \itemize{
#' \item The drifting semi-Markov kernel \eqn{q_{\frac{t}{n}}(u,v,l)},
#'   for every \eqn{t \in \{ 0, \dots, n \}} and \eqn{u \in E}, has its sums
#'   over \eqn{v} and \eqn{l}, equal to \eqn{1}:
#'
#'   \deqn{
#'   \sum_{v \in E}\sum_{l = 1}^{+\infty}q_{\frac{t}{n}}(u,v,l)
#'   = \sum_{v \in E}\sum_{l = 1}^{+\infty}A_{i}(t)\ q_{\frac{i}{d}}(u,v,l)
#'   = 1.}
#'
#' \item Therefore, we also get that for every \eqn{i \in \{0, \dots, d\}} and
#'   \eqn{u \in E}, the semi-Markov kernel \eqn{q_{\frac{i}{d}}(u,v,l)}
#'   has its sums over \eqn{v} and \eqn{l} equal to \eqn{1}:
#'   \deqn{\sum_{v \in E}\sum_{l = 1}^{+\infty}q_{\frac{i}{d}}(u,v,l)
#'   = 1.}
#'
#' \item Lastly, like in semi-Markov models, we do not allow sojourn times
#' equal to \eqn{0} or passing into the same state:
#' \deqn{q_{\frac{t}{n}}(u,v,0) = 0, \forall u,v \in E,}
#' \deqn{q_{\frac{t}{n}}(u,u,l) = 0, \forall u\in E,l\in\{1,\dots,+\infty\}.}
#' }
#'
#' \strong{Model specification restrictions}
#'
#' When we define a drifting semi-Markov model specification through the
#' functions \code{parametric_dsmm} or \code{nonparametric_dsmm},
#' the following restrictions need to be satisfied.
#'
#' \strong{\emph{Model 1}}
#'
#' The semi-Markov kernels are equal to
#' \eqn{q_{\frac{i}{d}}^{\ (1)}(u,v,l) =
#' p_{\frac{i}{d}}(u,v)f_{\frac{i}{d}}(u,v,l)}. Therefore,
#' \eqn{\forall u \in E} the sums
#' of \eqn{p_{\frac{i}{d}}(u,v)} over \eqn{v} and the sums of
#' \eqn{f_{\frac{i}{d}}(u,v,l)} over \eqn{l} must be
#' equal to 1:
#' \deqn{\sum_{v \in E} p_{\frac{i}{d}}(u,v) = 1,}
#' \deqn{\sum_{l = 1}^{+\infty }f_{\frac{i}{d}}(u,v,l) = 1.}
#'
#' \strong{\emph{Model 2}}
#'
#' The semi-Markov kernels are equal to \eqn{q_{\frac{i}{d}}^{\ (2)}(u,v,l) =
#'   p_{\frac{i}{d}}(u,v)f(u,v,l)}. Therefore, \eqn{\forall u \in E}
#'   the sums of \eqn{p_{\frac{i}{d}}(u,v)} over \eqn{v} and
#'   the sums of \eqn{f(u,v,l)} over \eqn{l} must be equal to 1:
#'   \deqn{\sum_{v \in E} p_{\frac{i}{d}}(u,v) = 1,}
#'   \deqn{\sum_{l = 1}^{+\infty }f(u,v,l) = 1.}
#'
#' \strong{\emph{Model 3}}
#'
#' The semi-Markov kernels are equal to \eqn{q_{\frac{i}{d}}^{\ (3)}(u,v,l) =
#'   p(u,v)f_{\frac{i}{d}}(u,v,l)}. Therefore,
#'   \eqn{\forall u \in E} the sums
#'   of \eqn{p(u,v)} over \eqn{v} and the sums of
#'   \eqn{f_{\frac{i}{d}}(u,v,l)} over \eqn{l} must be
#'   equal to 1:
#'   \deqn{\sum_{v \in E}p(u,v) = 1,}
#'   \deqn{\sum_{l = 1}^{+\infty }f_{\frac{i}{d}}(u,v,l) = 1.}
#'
#'
#' @section {Community Guidelines}:
#' For third parties wishing to contribute to the software, or to report issues
#' or problems about the software, they can do so directly through the
#' \href{https://github.com/Mavrogiannis-Ioannis/dsmmR}{development github page
#' of the package}.
#'
#'
#'
#' @section {Notes}:
#' Automated tests are in place in order to aid the user with any false input made
#' and, furthermore, to ensure that the functions used return the expected output.
#' Moreover, through strict automated tests, it is made possible for the user to
#' properly define their own \code{dsmm} objects and make use of them with the generic
#' functions of the package.
#'
#'
#' @seealso
#' For the estimation of a drifting semi-Markov model given a sequence:
#' \link{fit_dsmm}.
#'
#' For drifting semi-Markov model specifications:
#' \link{parametric_dsmm}, \link{nonparametric_dsmm}.
#'
#' For the simulation of sequences:
#' \link{simulate.dsmm}, \link{create_sequence}.
#'
#' For the retrieval of the drifting semi-Markov kernel through a
#' \code{dsmm} object: \link{get_kernel}.
#'
#' @keywords Drifting semi-Markov simulation estimation
#'
#' @references
#' Barbu, V. S., Limnios, N. (2008). Semi-Markov Chains and Hidden Semi-Markov
#' Models Toward Applications - Their Use in Reliability and DNA Analysis.
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#'
#' Vergne, N. (2008). Drifting Markov models with Polynomial Drift and
#' Applications to DNA Sequences. Statistical Applications in Genetics
#' Molecular Biology 7 (1).
#'
#' Barbu V. S., Vergne, N. (2019). Reliability and survival analysis for
#' drifting Markov models: modelling and estimation.
#' Methodology and Computing in Applied Probability, 21(4), 1407-1429.
#'
#' T. Nakagawa and S. Osaki. (1975). The discrete Weibull distribution.
#' IEEE Transactions on Reliability, R-24, 300-301.
#'
#' Sanger, F., Coulson, A. R., Hong, G. F., Hill, D. F., & Petersen, G. B.
#' (1982). Nucleotide sequence of bacteriophage \eqn{\lambda} DNA.
#' Journal of molecular biology, 162(4), 729-773.
#'
#' @author
#' \strong{Maintainer}: Ioannis Mavrogiannis \email{mavrogiannis.ioa@gmail.com}
#'
#' Authors:
#' \itemize{
#' \item Vlad Stefan Barbu
#' \item Ioannis Mavrogiannis
#' \item Nicolas Vergne
#' }
#'
#'
#'
#' @import stats
#' @import utils
#'
"_PACKAGE"
