#' @title dsmmR : Estimation and Simulation of Drifting Semi-Markov Models
#'
#' @description This package performs non-parametric estimation
#' and simulation based upon parametric and non-parametric model
#' specifications, for discrete-time Drifting Semi-Markov processes.
#' Three possible models are considered, regarding the number of transition
#' matrices and sojourn time distributions used for the computation of the
#' Semi-Markov kernels which characterize the Drifting Semi-Markov kernel.
#'
#' @details
#'
#' \strong{Introduction}
#'
#' The difference between the Markov models and the semi-Markov models concerns
#' the modeling of the sojourn time distributions.
#' The Markov Models (in discrete time) are modeled by a sojourn time
#' following the Geometric distribution. The Semi-Markov Models
#' are able to have a sojourn time distribution of arbitrary shape.
#' The further difference with a \emph{Drifting Semi-Markov Model},
#' is that we have \eqn{d + 1} (arbitrary) sojourn time distributions
#' and \eqn{d + 1} Transition Matrices (Model 1),
#' where \eqn{d} is defined as the polynomial degree.
#' Through them, we compute \eqn{d + 1} Semi-Markov kernels.
#' In this work, we also consider the possibility for obtaining these
#' Semi-Markov kernels with \eqn{d + 1} transition matrices and \eqn{1}
#' sojourn time distribution (Model 2) or \eqn{d + 1} sojourn time
#' distributions and \eqn{1} transition matrix (Model 3).
#'
#' @section Definition:
#'
#' Drifting Semi-Markov processes are particular non-homogeneous Markov
#' chains for which the Drifting Semi-Markov Kernel
#' \eqn{q_{\frac{t}{n}}(u,v,l)} is defined as
#' the probability of passing from the previous state \eqn{u} to the current
#' state \eqn{v} with a sojourn time of \eqn{l} at the instant \eqn{t}:
#' \deqn{q_{\frac{t}{n}}(u,v,l)=P(J_{t}=v,X_{t}=l|J_{t-1}=u),}
#' where \eqn{n} is the model size, defined as the length of the
#' embedded Markov chain and \eqn{X_{t}=S_{t}-S_{t-1}}
#' is the sojourn time at the instant \eqn{t}.
#'
#'
#' The Drifting Semi-Markov Kernel \eqn{q_{\frac{t}{n}}(u,v,l)} is a
#' linear combination of the product of \eqn{d + 1} Semi-Markov kernels:
#' \deqn{q_{\frac{i}{d}}(u,v,l)={p_{\frac{i}{d}}(u,v,l)}{f_{\frac{i}{d}}(u,v,l)},}
#' where for \eqn{i = 0,\dots,d} we have \eqn{d + 1} markov transition matrices
#' of the embedded Markov chain
#' \eqn{p_{\frac{i}{d}}(u,v,l)} and \eqn{d + 1} sojourn time distributions
#' \eqn{f_{\frac{i}{d}}(u,v,l)}.
#' This is the case where both \eqn{p} and \eqn{f}
#' are "drifting" between \eqn{d + 1} fixed points of the model, hence the
#' naming of "Drifting" Semi-Markov Models. We define the situation when
#' both \code{p} and \code{f} are drifting as Model 1.
#'
#' \strong{Model 1}
#'
#' Both \eqn{p} and \eqn{f} are drifting and the Drifting Semi-Markov Kernel
#' is described as:
#' \deqn{q_{\frac{t}{n}}(u,v,l) = \sum_{i = 0}^{d}A_{i}q_{\frac{i}{d}}(u,v,l) =
#' \sum_{i = 0}^{d}A_{i}{p_{\frac{i}{d}}(u,v,l)}{f_{\frac{i}{d}}(u,v,l)},} where
#' \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree \eqn{d},
#' that satisfy the conditions:
#' \eqn{\forall t, \sum_{i=0}^{d}A_{i}(t) = 1} and \eqn{A_i(\frac{nj}{d})=
#' \mathbf{1}_{\{i=j\}}}, with \eqn{\mathbf{1}_{\{i=j\}}} the function
#' that is \eqn{1} when \eqn{i = j} and \eqn{0} otherwise.
#'
#' \strong{Model 2}
#'
#' In this case, \eqn{p} is drifting and \eqn{f} is \strong{not drifting}.
#' Therefore, the Drifting Semi-Markov Kernel is now described as:
#' \deqn{q_{\frac{t}{n}}(u,v,l) = \sum_{i = 0}^{d}A_{i}q_{\frac{i}{d}}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}{p_{\frac{i}{d}}(u,v,l)}{f(u,v,l)},}
#' where \eqn{f(u,v,l)} remains constant for all instances \eqn{t}
#' alongside the embedded Markov chain with length \eqn{n}.
#'
#'
#' \strong{Model 3}
#'
#' In this case, \eqn{f} is drifting and \eqn{p} is \strong{not drifting}.
#' Therefore, the Markov Kernel is now described as:
#' \deqn{q_{\frac{t}{n}}(u,v,l) = \sum_{i = 0}^{d}A_{i}q_{\frac{i}{d}}(u,v,l)
#' = \sum_{i = 0}^{d}A_{i}{p(u,v,l)}{f_{\frac{i}{d}}(u,v,l)},}
#' where \eqn{p(u,v)} remains constant for all instances \eqn{t}
#' alongside the embedded Markov chain with length \eqn{n}.
#'
#'
#' @section Parametric and Non-parametric model specifications :
#'
#' In this package, we can define parametric and non-parametric Drifting
#' Semi-Markov models.
#' For the \emph{parametric} case, several discrete distributions are
#' considered for the modeling of the sojourn times:
#' Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial.
#' They can be specified through the function
#' \code{parametric_dsmm} that returns an object with the
#' S3 class (\code{dsmm}, \code{dsmm_parametric}).
#'
#' The \emph{non-parametric} model specification concerns the sojourn
#' time distributions when \strong{no assumptions} are done about the
#' shape of the distributions. This is possible through the function
#' \code{nonparametric_dsmm()}, that returns  an object of class
#' (\code{dsmm}, \code{dsmm_nonparametric}).
#'
#' It is also possible to proceed with a non-parametric estimation
#' for a model on an existing sequence, through the function
#' \code{fit_dsmm()}, which returns an object with the S3 class
#' (\code{dsmm}, \code{dsmm_fit}).
#'
#'
#' Based on an \code{dsmm} object, it is possible to use the following methods:
#' \itemize{
#'   \item Simulate a sequence through the function \code{simulate.dsmm()};
#'   \item Get the Drifting Semi-Markov Kernel, \eqn{q_{\frac{t}{n}}(u,v,l)}
#'    by using the generic function \code{get_kernel()}.
#' }
#'
#' These are further specified when necessary, depending on the classes
#' \code{dsmm_fit}, \code{dsmm_parametric} and \code{dsmm_nonparametric}.
#' In summary, the \code{dsmm} class acts like a wrapper class
#' for Drifting Semi-Markov model specifications, while the classes
#' \code{dsmm_parametric}, \code{dsmm_nonparametric} and \code{dsmm_fit}
#' are exclusive to the functions that create the corresponding models,
#' and inherit methods from the \code{dsmm} class.
#'
#' @section Restrictions :
#'
#' The following restrictions must be satisfied for every defined
#' Drifting Semi-Markov Model:
#' \itemize{
#' \item The Drifting Semi-Markov kernel \eqn{q_{\frac{t}{n}}(u,v,l)},
#'   for every \eqn{t \in \{ 0, \dots, n \}} and \eqn{u \in E}, has its sums
#'   over \eqn{v} and over \eqn{l} equal to \eqn{1}:
#'   \deqn{{\sum_{v \in E}{\sum_{l = 1}^{k_{max}}q_{\frac{t}{n}}(u,v,l)}} = 1.}
#'
#' \item Furthermore, since the Drifting Semi-Markov kernel
#'   \eqn{q_{\frac{t}{n}}(u,v,l) =
#'     \sum_{i = 0}^{d}\sum_{t=0}^{n}A_{i}(t)q_{\frac{t}{n}}(u,v,l)},
#'   we also get that for every \eqn{i \in \{0, \dots, d\}} and
#'   \eqn{u \in E}, the kernel \eqn{q_{\frac{i}{d}}(u,v,l)} has its sums
#'   over \eqn{v} and over \eqn{l} equal to \eqn{1}:
#'   \deqn{\sum_{v \in E}\sum_{l = 1}^{k_{max}}q_{\frac{t}{n}}(u,v,l) = 1.}
#'
#' \item Lastly, we do not allow sojourn times equal to \eqn{0} or passing into
#' the same state:
#' \deqn{q_{\frac{t}{n}}(u,v,0) = 0,}
#' \deqn{q_{\frac{t}{n}}(u,u,l) = 0.}
#' }
#'
#'
#'
#' Specifying the models, the following are necessary conditions:
#'
#' \strong{Model 1}
#'
#' The Semi-Markov kernels \eqn{q_{\frac{i}{d}}(u,v,l) =
#'   p_{\frac{i}{d}}(u,v)f_{\frac{i}{d}}(u,v,l)}. Therefore,
#'  \eqn{\forall t \in \{ 0, \dots, n \}} and \eqn{\forall u \in E}, the sums
#'   of \eqn{p_{\frac{i}{d}}(u,v)} over \eqn{v\in E}, and the sums of
#'   \eqn{f_{\frac{i}{d}}(u,v,l)} over \eqn{l\in \{0,\dots,k_{max} \}}, must be
#'   equal to 1:
#'   \deqn{\sum_{v \in E} p_{\frac{i}{d}}(u,v) = 1,}
#'   \deqn{\sum_{l = 1}^{k_{max}}f_{\frac{i}{d}}(u,v,l) = 1.}
#'
#' \strong{Model 2}
#'
#' The Semi-Markov kernels \eqn{q_{\frac{i}{d}}(u,v,l) =
#'   p_{\frac{i}{d}}(u,v)f(u,v,l)}. Therefore,
#'   \eqn{\forall t \in \{ 0, \dots, n \}} and \eqn{\forall u \in E}, the sums
#'   of \eqn{p_{\frac{i}{d}}(u,v)} over \eqn{v\in E}, and the sums of
#'   \eqn{f(u,v,l)} over \eqn{l\in\{0,\dots,k_{max}\}}, must be
#'   equal to 1:
#'   \deqn{\sum_{v \in E} p_{\frac{i}{d}}(u,v) = 1,}
#'   \deqn{\sum_{l = 1}^{k_{max}}f(u,v,l) = 1.}
#'
#' \strong{Model 3}
#'
#' The Semi-Markov kernels \eqn{q_{\frac{i}{d}}(u,v,l) =
#'   p(u,v)f_{\frac{i}{d}}(u,v,l)}. Therefore,
#'   \eqn{\forall t \in \{ 0, \dots, n \}} and \eqn{\forall u \in E}, the sums
#'   of \eqn{p(u,v)} over \eqn{v\in E}, and the sums of
#'   \eqn{f_{\frac{i}{d}}(u,v,l)} over \eqn{l\in\{0,\dots,k_{max}\}}, must be
#'   equal to 1:
#'   \deqn{\sum_{v \in E}p(u,v) = 1,}
#'   \deqn{\sum_{l = 1}^{k_{max}}f_{\frac{i}{d}}(u,v,l) = 1.}
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
