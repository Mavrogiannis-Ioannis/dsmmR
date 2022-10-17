#' @title dsmmR : Simulation and Estimation of Drifting Semi-Markov Models
#'
#' @description This package performs parametric and non-parametric estimation
#' and simulation for multi-state discrete-time Drifting Semi-Markov processes,
#' for 3 possible Model specifications.
#'
#' @details
#'
#' \strong{Introduction}
#'
#' The difference between the Markov model and the semi-Markov model concerns
#' the modelling of the sojourn time.
#' The Markov Models are modeled by a (discrete time) Geometric distribution.
#' However, for a Semi-Markov Model, the sojourn time can be \emph{arbitrarily}
#' distributed. The further difference with a \emph{Drifting Semi-Markov Model},
#' is that we have \eqn{d + 1} Sojourn Time Distributions that are
#' modeled after arbitrary discrete distributions, where \eqn{d} is the
#' polynomial degree of the drift.
#'
#' @section Definition:
#'
#' Drifting Semi-Markov processes are particular non-homogeneous Markov
#' chains, for which the Drifting Semi-Markov Kernel
#' \eqn{q_{\frac{t}{n}}(u,v,l)} is defined as
#' the probability of passing from the previous state \eqn{u} to the current
#' state \eqn{v}, when the sojourn time spent at \eqn{u} is \eqn{l},
#' at the instant \eqn{t}:
#' \deqn{q_{\frac{t}{n}}(u,v,l)=P(J_{t}=v,X_{t}=l|J_{t-1}=u),}
#' where \eqn{X_{t}=S_{t}-S_{t-1}}, is the sojourn time at the instant \eqn{t}.
# Also, the probability that the current state \eqn{v}
# will be equal to the previous state \eqn{i} is 0, for all \eqn{l}
#'
#'
#' The Drifting Semi-Markov Kernel \eqn{q_{\frac{t}{n}}(u,v,l)} is a
#' linear combination of the product of \eqn{d + 1} kernels:
#' \deqn{q_{\frac{i}{d}}(u,v,l)={p_{\frac{i}{d}}(u,v,l)}{f_{\frac{i}{d}}(u,v,l)},}
#' where for \eqn{i = 0,\dots,d} we have \eqn{d + 1} Markov Transition Matrices
#' of the embedded Markov chain,
#' \eqn{p_{\frac{i}{d}}(u,v,l)}, and \eqn{d + 1} Sojourn Time Distributions
#' \eqn{f_{\frac{i}{d}}(u,v,l)}, for \eqn{d} the polynomial degree.
#' This is the case where both \eqn{p} and \eqn{f}
#' are "drifting" between \eqn{d + 1} fixed points of the model, hence the
#' naming of "Drifting" Semi-Markov Models. When both \code{p} and \code{f} are
#' drifting we have the first model, Model 1.
#'
#' \strong{Model 1}
#'
#' When both \eqn{p} and \eqn{f} are drifting, the Drifting Semi-Markov Kernel
#' is described as:
#' \deqn{q_{\frac{t}{n}}(u,v,l) = \sum_{i = 0}^{d}A_{i}q_{\frac{i}{d}}(u,v,l) =
#' \sum_{i = 0}^{d}A_{i}{p_{\frac{i}{d}}(u,v,l)}{f_{\frac{i}{d}}(u,v,l)},} where
#' \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree \eqn{d},
#' that satisfy certain conditions, most notably that
#' \eqn{\forall t, \sum_{i=0}^{d}A_{i}(t) = 1}.
#' Two more models are considered in this work:
#'
#' \strong{Model 2}
#'
#' In this case, \eqn{p} is drifting and \eqn{f} is \strong{not drifting}.
#' Therefore, the Drifting Semi-Markov Kernel is now described as:
#' \deqn{q_{\frac{t}{n}}(u,v,l) = \sum_{i = 0}^{d}A_{i}q_{\frac{i}{d}}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}{p_{\frac{i}{d}}(u,v,l)}{f(u,v,l)},}
#' where \eqn{f(u,v,l)} is \strong{not drifting} alongside the point \eqn{t}
#' of the model with size (or length) \eqn{n}.
#'
#'
#' \strong{Model 3}
#'
#' In this case, \eqn{f} is drifting and \eqn{p} is \strong{not drifting}.
#' Therefore, the Markov Kernel is now described as:
#' \deqn{q_{\frac{t}{n}}(u,v,l) = \sum_{i = 0}^{d}A_{i}q_{\frac{i}{d}}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}{p(u,v,l)}{f_{\frac{i}{d}}(u,v,l)},}
#' where \eqn{p(u,v)} is \strong{not drifting} alongside the point \eqn{t}
#' of the model with size (or length) \eqn{n}.
#'
#'
#' @section Parametric and Non-parametric modelling :
#'
#' For the \emph{parametric estimation} of a Drifting Semi-Markov Model,
#' several discrete distributions are considered for the
#' sojourn times:
#' Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial.
#' The \emph{non-parametric estimation} concerns the sojourn time
#' distributions when \strong{no assumptions} are done about the shape of
#' distributions.
#' They can be specified through the function
#' \code{parametric_dsmm} that returns an object with the
#' S3 class (\code{dsmm}, \code{dsmm_parametric}).
#'
#' Regarding the non-parametric specification, it is possible to define a model
#' through the function \code{nonparametric_dsmm()} or to fit a model on an
#' existing sequence, through the function \code{fit_dsmm()}.
#' These functions return objects with the S3 classes
#' (\code{dsmm}, \code{dsmm_nonparametric}) and (\code{dsmm}, \code{dsmm_fit}),
#' respectively.
#'
#' This means that the \code{dsmm} class acts like a wrapper class
#' (or a parent class) for Drifting Semi-Markov model specifications,
#' while the classes
#' \code{dsmm_parametric}, \code{dsmm_nonparametric} and \code{dsmm_fit}
#' are exclusive to the functions that create the corresponding models,
#' and inherit methods from the \code{dsmm} class.
#'
#' Based on a model specification through the aforementioned functions,
#' it is possible to use the following methods:
#' \itemize{
#'   \item Simulate a sequence through the function \code{simulate.dsmm()};
#'   \item Get the Drifting Semi-Markov Kernel, \eqn{q_{\frac{t}{n}}(u,v,l)}
#'    by using the generic function \code{get_kernel()}.
#' }
#'
#'
#' @section Restrictions :
#'
#' The following restrictions must be satisfied for every defined
#' Drifting Semi-Markov Model:
#' \itemize{
#' \item The Drifting Semi-Markov kernel \eqn{q_{\frac{t}{n}}(u,v,l)},
#'   for every \eqn{t \in \{ 0, \dots, n \}} and \eqn{u \in E}, has its sums
#'   over \eqn{v} and over \eqn{l} equal to \eqn{1}:
#'   \deqn{{\sum_{v \in E}{\sum_{l = 1}^{k_{max}}q_{\frac{t}{n}}(u,v,l)}} = 1}
#'
#' \item Furthermore, since the Drifting Semi-Markov kernel
#'   \eqn{q_{\frac{t}{n}}(u,v,l) =
#'     \sum_{i = 0}^{d}\sum_{t=0}^{n}A_{i}(t)q_{\frac{t}{n}}(u,v,l)},
#'   we also get that for every \eqn{i \in \{0, \dots, d\}} and
#'   \eqn{u \in E}, the kernel \eqn{q_{\frac{i}{d}}(u,v,l)} has its sums
#'   over \eqn{v} and over \eqn{l} equal to \eqn{1}:
#'   \deqn{\sum_{v \in E}\sum_{l = 1}^{k_{max}}q_{\frac{t}{n}}(u,v,l) = 1}
#' }
#'
#' Specifying the models, the following are necessary conditions:
#'
#' \strong{Model 1}
#'
#' The kernel \eqn{q_{\frac{i}{d}}(u,v,l) =
#'   p_{\frac{i}{d}}(u,v)f_{\frac{i}{d}}(u,v,l)}. Therefore,
#'  \eqn{\forall t \in \{ 0, \dots, n \}} and \eqn{\forall u \in E}, the sums
#'   of \eqn{p_{\frac{i}{d}}(u,v)} over \eqn{v\in E}, and the sums of
#'   \eqn{f_{\frac{i}{d}}(u,v,l)} over \eqn{l\in \{0,\dots,k_{max} \}}, must be
#'   equal to 1:
#'   \deqn{\sum_{v \in E} p_{\frac{i}{d}}(u,v) = 1}
#'   \deqn{\sum_{l = 1}^{k_{max}}f_{\frac{i}{d}}(u,v,l) = 1}
#'
#' \strong{Model 2}
#'
#' The kernel \eqn{q_{\frac{i}{d}}(u,v,l) =
#'   p_{\frac{i}{d}}(u,v)f(u,v,l)}. Therefore,
#'   \eqn{\forall t \in \{ 0, \dots, n \}} and \eqn{\forall u \in E}, the sums
#'   of \eqn{p_{\frac{i}{d}}(u,v)} over \eqn{v\in E}, and the sums of
#'   \eqn{f(u,v,l)} over \eqn{l\in\{0,\dots,k_{max}\}}, must be
#'   equal to 1:
#'   \deqn{\sum_{v \in E} p_{\frac{i}{d}}(u,v) = 1}
#'   \deqn{\sum_{l = 1}^{k_{max}}f(u,v,l) = 1}
#'
#' \strong{Model 3}
#'
#' The kernel \eqn{q_{\frac{i}{d}}(u,v,l) =
#'   p(u,v)f_{\frac{i}{d}}(u,v,l)}. Therefore,
#'   \eqn{\forall t \in \{ 0, \dots, n \}} and \eqn{\forall u \in E}, the sums
#'   of \eqn{p(u,v)} over \eqn{v\in E}, and the sums of
#'   \eqn{f_{\frac{i}{d}}(u,v,l)} over \eqn{l\in\{0,\dots,k_{max}\}}, must be
#'   equal to 1:
#'   \deqn{\sum_{v \in E}p(u,v) = 1}
#'   \deqn{\sum_{l = 1}^{k_{max}}f_{\frac{i}{d}}(u,v,l) = 1}
#'
#' @seealso \link{fit_dsmm}, \link{parametric_dsmm}, \link{nonparametric_dsmm},
#' \link{simulate}, \link{get_kernel}
#'
#' @keywords Drifting Semi-Markov simulation estimation
#'
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov
#' Models Toward Applications - Their Use in Reliability and DNA Analysis.
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#'
#' N. Vergne. (2008). Drifting Markov Models with Polynomial Drift and
#' Applications to DNA Sequences. Statistical Applications in Genetics
#' Molecular Biology 7 (1).
#'
#' T. Nakagawa and S. Osaki. (1975). The discrete Weibull distribution.
#' IEEE Transactions on Reliability, R-24, 300-301.
#'
#' @import stats
#' @import utils
#'
"_PACKAGE"
