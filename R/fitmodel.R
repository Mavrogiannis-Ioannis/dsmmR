# '''
#    1. This file concerns itself with the creation and definition of the
#    fitted drifting semi-Markov model on a random sequence, using `fit_dsmm()`.
#    It returns two possible classes: `dsmm_fit_parametric` or
#    `dsmm_fit_nonparametric`. These are child classes of `dsmm`.
#    They inherit the methods:
#    `simulate`, `get_kernel`.
#    It has its own methods for:
#    `print`, `check_attributes`.
# '''


#' @title Estimation of a Drifting semi-Markov chain
#' @aliases dsmm_fit dsmm_fit_nonparametric dsmm_fit_parametric
#' @description Least Squares Estimation
#' of a Drifting semi-Markov chain, given one sequence of states.
#' This estimation can be parametric or non-parametric and is available for
#' the three types of Drifting semi-Markov models.
#'
#' @param sequence Character vector that represents a sequence of states in
#'     \eqn{E}. States must be characters with \eqn{length \geq 1}.
#' @param states Character vector that represents the state space
#'     \eqn{E} of choice, with length equal to \eqn{|E| = s}.
#' @param degree Positive integer that represents the polynomial degree
#'     \eqn{d} for the Drifting semi-Markov Model.
#' @param f_is_drifting Logical. Specifies if \eqn{f} is drifting or not.
#' @param p_is_drifting Logical. Specifies if \eqn{p} is drifting or not.
#' @param initial_dist Optional. Character that describes the method to
#'     estimate the initial distribution.
#'     \itemize{
#'     \item
#'           \code{"unif"} : the initial probability of each state is
#'           equal to \eqn{1/s}.
#'     \item
#'           \code{"freq"} : the initial distribution of each state is
#'           equal to the frequencies of each state in the sequence.
#'     }
#' @param estimation Optional. Character vector. Specifies whether the
#' .    estimation will be nonparametric or parametric.
#'      Possible values are "nonparametric" or "parametric".
#' @param f_dist Optional. Can be defined in two ways:
#'     \itemize{
#'     \item
#'          If \code{estimation} is equal to "nonparametric", it is equal to
#'          \code{NULL}.
#'     \item
#'         if \code{estimation} is equal to "parametric", it is a
#'         character array that specifies the distributions of the sojourn
#'         times, for every state transition.
#'         The list of possible values is:
#'         \code{["unif", "geom", "pois", "dweibull", "nbinom", NA]}.
#'         It can be defined in two ways:
#'         \itemize{
#'             \item
#'                 If \eqn{f} \strong{is not} drifting, it has dimensions of
#'                 \eqn{(s \times s)}.
#'             \item
#'                 If \eqn{f} \strong{is} drifting, it has dimensions of
#'                 \eqn{(s \times s \times d+1)}
#'                 (see more in \emph{Details, Defined Arguments}.)
#'         }
#'         It is defined similarly to the attribute \code{f_dist}
#'         in \link{dsmm_parametric}.
#'     }
# @param numerical_est Optional. Logical. Specifies if numerical estimation
#   under constraint should be used instead of LSE.
#    \emph{Currently not supported}. Default value is \code{FALSE}.
#'
#' @details This function estimates a Drifting semi-Markov Model in the
#'     parametric and non-parametric case.
#'     The parametric estimation can be achieved by fitting the non-parametric
#'     estimations of the sojourn time distributions, with regards to the
#'     distributions chosen by the user in the attribute \code{f_dist}.
#'     It is possible for three different models to be estimated
#'     (more about them in \link{dsmmR}).
#'     A normalization technique is used in order to correct estimation errors
#'     from small sequences.
#'
#' \strong{Non-parametric Estimation}
#'
#' We try to minimize the value.(maybe not?)
#'
#' \strong{\emph{Model 1}}
#'
#' When the transition matrix of the embedded Markov chain \eqn{p} and
#' the conditional sojourn time distribution \eqn{f} are both drifting,
#' the Drifting semi-Markov kernel can be estimated as:
#' \deqn{\hat{q}_{\frac{t}{n}}^{(1)}(u,v,l) =
#'       \sum_{i = 0}^{d}A_{i}(t)\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l),}
#' \eqn{\forall t \in \{0,\dots,n\}, \forall u,v\in E,
#' \forall l \in \{0,\dots, k_{max} \} }, where \eqn{k_{max}} is the maximum
#' sojourn time that was observed in the sequence and
#' \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree
#' \eqn{d} (see \link{dsmmR}).
#'
#' The semi-Markov kernels
#' \eqn{\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l), i = 0, \dots, d},
#' are estimated through Least Squares Estimation and are obtained
#' after solving the system
#'
#' \deqn{MJ = P \iff \\
#' \left(\begin{array}{cc}
#'       \sum_{t=1}^{n}1_{u}(t)A_{0}(t)A_{0}(t) & \dots &
#'       \sum_{t=1}^{n}1_{u}(t)A_{0}(t)A_{d}(t)\\
#'       \vdots & \ddots & \vdots \\
#'       \sum_{t=1}^{n}1_{u}(t)A_{d}(t)A_{0}(t) & \dots &
#'       \sum_{t=1}^{n}1_{u}(t)A_{d}(t)A_{d}(t)
#'  \end{array}\right)
#'  \left(\begin{array}{cc}
#'       \hat{q}_{0}^{(1)}(u,v,l) \\
#'       \vdots \\
#'       \hat{q}_{\frac{i}{d}}^{(1)}(u,v,l) \\
#'       \vdots \\
#'       \hat{q}_{1}^{(1)}(u,v,l)
#'  \end{array} \right)
#'       =
#'  \left(\begin{array}{cc}
#'       \sum_{t=1}^{n}1_{uvl}(t)A_{0}(t) \\
#'       \vdots \\
#'       \sum_{t=1}^{n}1_{uvl}(t)A_{i}(t) \\
#'       \vdots \\
#'       \sum_{t=1}^{n}1_{uvl}(t)A_{d}(t)
#'  \end{array} \right)}
#' where
#' \itemize{
#' \item
#' \eqn{M = (M_{ij})_{i,j \in \{0, \dots, d\} } =
#' (\sum_{t=1}^{n}1_{u}(t)A_{i}(t)A_{j}(t))_{
#' i,j \in \{0, \dots, d\}}}
#'
#' \item
#' \eqn{J = (J_i)_{i \in \{0, \dots, d\} } =
#' \left(\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)\right)_{
#' i \in \{0, \dots, d\}}
#' }
#'
#' \item
#' \eqn{P=(P_i)_{i\in \{0, \dots, d\} }=
#' (\sum_{t=1}^{n}1_{uvl}(t)A_{i}(t))_{
#' i \in \{0, \dots, d\}}
#' }
#' }
#' where we use the indicator functions:
#' \deqn{1_{u}(t) =
#' \begin{cases}
#'   1  & \text{if at } t \text{ the previous state is } u, \\
#'   0  & \text{otherwise}.
#' \end{cases}
#' }
#' and
#' \deqn{1_{uvl}(t) =
#' \begin{cases}
#'  1  & \text{if at } t \text{ the previous state is } u
#'  \text{, with sojourn
#'  time } l \text{ and next state } v, \\
#'  0  & \text{otherwise}.
#' \end{cases}
#' }
#'
#' In order to obtain the estimations of \eqn{\hat{p}_{\frac{i}{d}}(u,v)}
#' and \eqn{\hat{f}_{\frac{i}{d}}(u,v,l)}, we use the following estimations:
#'
#' \deqn{\hat{p}_{\frac{i}{d}}(u,v) =
#'     \sum_{l = 0}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l),}
#' \deqn{\hat{f}_{\frac{i}{d}}(u,v,l) =
#'     \frac{\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}{
#'          \sum_{l = 0}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}.}
#'
#'
#' \strong{\emph{Model 2}}
#'
#' In this case, \eqn{p} is drifting and \eqn{f} is not drifting. Therefore,
#' the estimated Drifting semi-Markov kernel will be given by:
#' \deqn{\hat{q}_{\frac{t}{n}}^{(2)}(u,v,l) =
#' \sum_{i=0}^{d}\hat{q}_{\frac{i}{d}}^{(2)}(u,v,l),}
#' \eqn{\forall t \in \{0,\dots,n\}, \forall u,v\in E,
#' \forall l\in \{0,\dots, k_{max} \}},
#' where \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree
#' \eqn{d} (see \link{dsmmR}). Since \eqn{p} is drifting,
#' we define the estimation of \eqn{p} the same way as we did in Model 1.
#' Therefore,
#' \eqn{\forall u,v \in E, \forall l \in \{0,\dots, k_{max} \}},
#' we have the following estimations:
#'
#' \deqn{\hat{p}_{\frac{i}{d}}(u,v) =
#'     \sum_{l = 0}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l),}
#' \deqn{\hat{f}(u,v,l) =
#'     \frac{\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}{
#'        \sum_{i = 0}^{d}\sum_{l=0}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)},}
#'
#' Thus, the \emph{estimated} semi-Markov kernels for Model 2,
#' \eqn{\hat{q}_{\frac{i}{d}}^{(2)}(u,v,l) =
#' \hat{p}_{\frac{i}{d}}(u,v)\hat{f}(u,v,l)}, can be written with
#' regards to the \emph{estimated} semi-Markov kernels of Model 1,
#' \eqn{\hat{q}_{\frac{i}{d}}}, as in the following:
#'
#' \deqn{\hat{q}_{\frac{i}{d}}^{(2)}(u,v,l) = \frac{
#' \sum_{l=1}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)
#' \sum_{i = 0}^{d}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}{
#' \sum_{i = 0}^{d}\sum_{l=1}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}.}
#'
#'
#' \strong{\emph{Model 3}}
#'
#' In this case, \eqn{f} is drifting and \eqn{p} is not drifting. Therefore,
#' the estimated Drifting semi-Markov kernel will be given by:
#' \deqn{\hat{q}_{\frac{t}{n}}^{(3)}(u,v,l) =
#' \sum_{i=0}^{d}\hat{q}_{\frac{i}{d}}^{(3)}(u,v,l),}
#' \eqn{\forall t \in \{0,\dots,n\}, \forall u,v\in E,
#' \forall l\in \{0,\dots, k_{max} \}},
#' where \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree
#' \eqn{d} (see \link{dsmmR}). Since \eqn{f} is drifting,
#' we define the estimation of \eqn{f} the same way as we did in Model 1.
#' Therefore,
#' \eqn{\forall u,v \in E, \forall l \in \{0,\dots, k_{max} \}},
#' we have the following estimations:
#'
#' \deqn{\hat{p}(u,v) =
#' \frac{\sum_{i=0}^{d}\sum_{l=0}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}{d+1}}
#' \deqn{\hat{f}_{\frac{i}{d}}(u,v,l) =
#'     \frac{\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}{
#'          \sum_{l = 0}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}.}
#'
#' Thus, the \emph{estimated} semi-Markov kernels for Model 3,
#' \eqn{\hat{q}_{\frac{i}{d}}^{(3)}(u,v,l) =
#' \hat{p}(u,v)\hat{f}_{\frac{i}{d}}(u,v,l)}, can be written with
#' regards to the \emph{estimated} semi-Markov kernels of Model 1,
#' \eqn{\hat{q}_{\frac{i}{d}}}, as in the following:
#'
#' \deqn{\hat{q}_{\frac{i}{d}}^{(3)}(u,v,l) = \frac{
#' \hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)
#' \sum_{i=0}^{d}\sum_{l=0}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}
#' {(d+1)\sum_{l=1}^{k_{max}}\hat{q}_{\frac{i}{d}}^{(1)}(u,v,l)}.}
#'
#'
#' @return Returns an object of S3 class \code{(dsmm_fit_nonparametric, dsmm)} or
#' \code{(dsmm_fit_parametric, dsmm)}.
#' It has the following attributes:
#' \itemize{
#' \item \code{dist} : List that contains the \emph{p} and \emph{f}
#' estimated distributions from the sequence, given in the argument
#' \code{sequence};
#' \item \code{seq} : Character vector that contains the
#' \strong{state jumps of the original sequence}. It is this attribute of the
#' object that describes the length of the model \eqn{n}. Last state is also
#' included, for a total length of \eqn{n+1}, but it should not used;
#' \item \code{soj_times} : Numerical vector that contains the sojourn times
#' spent for each state in \code{seq} before the jump to the next state;
#' Last state is also
#' included, for a total length of \eqn{n+1}, but it should not used;
#' \item \code{k_max} : Numerical value that contains the maximum sojourn
#' time, meaning the maximum value in \code{soj_times};
#' \item \code{initial_dist} : Numerical vector that contains an estimation
#' for the initial distribution of the realized states in \code{sequence};
#' \item \code{model_size} : Integer value that contains the length of the
#' model This is equal to \code{length(seq) - 1}, for \code{seq} as defined
#' above. \item \code{states} : Character vector that contains the realized
#' states given in the argument \code{sequence};
#' \item \code{s} : Integer that contains the length of the state
#' space \eqn{E} given in the attribute \code{states}.
#' \item \code{degree} : Integer that contains the polynomial degree
#' \eqn{d} considered for the drifting of the model.
#' \item \code{f_is_drifting} : Logical. Passing down from the arguments.
#' \item \code{p_is_drifting} : Logical. Passing down from the arguments.
# \item \code{numerical_est} : Logical. Passing down from the arguments.
#' \item \code{Model} : Character vector. Possible values:
#' \code{["Model 1", "Model 2", "Model 3"]}, corresponding to whether
#' \eqn{p,f} are drifting or not.
#' \item \code{estimation} : Character vector. Specifies whether parametric or
#' nonparametric estimation was used.
#' \item \code{A_i} : Numerical Matrix. Used for the methods defined for the
#' object. Is not printed when viewing the object.
#' \item \code{Ji} : Numerical Array. Used for the methods defined for the
#' object. Is not printed when viewing the object.
#' }
#'
#' @seealso
#' For the theoretical background of Drifting semi-Markov Models: \link{dsmmR}.
#'
#' For sequence simulation: \link{simulate.dsmm} and \link{create_sequence}.
#'
#' For Drifting semi-Markov model specification:
#' \link{parametric_dsmm}, \link{nonparametric_dsmm}
#'
#' For the retrieval of the Drifting semi-Markov kernel:
#' \link{get_kernel}.
#'
#' @references
#' V. S. Barbu, N. Limnios. (2008). semi-Markov Chains and Hidden semi-Markov
#' Models Toward Applications - Their Use in Reliability and DNA Analysis.
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#'
#'
#'
#' @export
#'
#' @examples
#' # Create a random sequence
#' sequence <- create_sequence("DNA", len = 2000, seed = 1)
#' # Alternatively, we could use ``
#' # Alternatively, we could use `create_sequence("DNA")`
#' # to retrieve a randomly generated sequence.
#' # > data("lambda", package = "dsmmR")
#' # > sequence <- c(lambda)
#' states <- sort(unique(sequence))
#' degree <- 3
#'
#' # ===========================================================================
#' # Nonparametric Estimation.
#' # Fitting the `lambda` genome under distributions of unknown shape.
#' # ===========================================================================
#'
#' # ---------------------------------------------------------------------------
#' # Both p and f are drifting - Model 1.
#' obj_model_1 <- fit_dsmm(sequence = sequence,
#'                        states = states,
#'                        degree = degree,
#'                        f_is_drifting = TRUE,
#'                        p_is_drifting = TRUE,
#'                        initial_dist = "freq",
#'                        estimation = "nonparametric", # default value
#'                        f_dist = NULL # default value
#'                        )
#' cat(paste0(
#'     "We fitted a sequence with ", obj_model_1$Model, ",\n",
#'     "model size: n = ", obj_model_1$model_size, ",\n",
#'     "length of state space: s = ", obj_model_1$s, ",\n",
#'     "maximum sojourn time: k_max = ", obj_model_1$k_max, " and\n",
#'     "polynomial (drifting) Degree: d = ", obj_model_1$degree, ".\n"
#' ))
#' # Get the drifting p and f arrays.
#' p_drift <- obj_model_1$dist$p_drift
#' f_drift <- obj_model_1$dist$f_drift
#' cat(paste0(
#'     "Dimension of p_drift: (s, s, d + 1) = (",
#'     paste(dim(p_drift), collapse = ", "), ").\n",
#'     "Dimension of f_drift: (s, s, k_max, d + 1) = (",
#'     paste(dim(f_drift), collapse = ", "), ").\n"
#' ))
#'
#' # ---------------------------------------------------------------------------
#' # Fitting the sequence when p is drifting and f is not drifting - Model 2.
#' obj_model_2 <- fit_dsmm(sequence = sequence,
#'                        states = states,
#'                        degree = degree,
#'                        initial_dist = "unif",
#'                        f_is_drifting = FALSE,
#'                        p_is_drifting = TRUE)
#' cat(paste0("We fitted a sequence with ", obj_model_2$Model, ".\n"))
#' # Get the drifting p and non-drifting f arrays.
#' p_drift_2 <- obj_model_2$dist$p_drift
#' f_notdrift <- obj_model_2$dist$f_notdrift
#' all.equal.numeric(p_drift, p_drift_2) # p is the same as in Model 1.
#' cat(paste0(
#'     "Dimension of f_notdrift: (s, s, k_max) = (",
#'      paste(dim(f_notdrift), collapse = ", "), ").\n"
#' ))
#'
#' # ---------------------------------------------------------------------------
#' # Fitting the sequence when f is drifting and p is not drifting - Model 3.
#' # ---------------------------------------------------------------------------
#'
#' obj_model_3 <- fit_dsmm(sequence = sequence,
#'                        states = states,
#'                        degree = degree,
#'                        f_is_drifting = TRUE,
#'                        p_is_drifting = FALSE)
#' cat(paste0("We fitted a sequence with ", obj_model_3$Model, ".\n"))
#' # Get the drifting f and non-drifting p arrays.
#' p_notdrift <- obj_model_3$dist$p_notdrift
#' f_drift_3 <- obj_model_3$dist$f_drift
#' all.equal.numeric(f_drift, f_drift_3) # f is the same as in Model 1.
#' cat(paste0("Dimension of f_notdrift: (s, s) = (",
#'     paste(dim(p_notdrift), collapse = ", "), ").\n"))
#'
#'
#' # ===========================================================================
#' # Parametric Estimation
#' # Fitting the `lambda` genome under distributions of known shape.
#' # ===========================================================================
#' ### Comments
#' ### 1.  For the parametric estimation it is recommended to use a common set
#' ###     of distributions while only the parameters (of the sojourn times)
#' ###     are drifting. This results in (generally) higher accuracy.
#' ### 2.  This process is similar to that used in `dsmm_parametric()`.
#'
#' s <- length(states)
#' # Getting the distributions for the states.
#' # Rows correspond to previous state `u`.
#' # Columns correspond to next state `v`.
#' f_dist_1 <- matrix(c(NA,         "unif",     "dweibull", "nbinom",
#'                      "pois",      NA,        "pois",     "dweibull",
#'                      "geom",     "pois",      NA,        "geom",
#'                      "dweibull", 'geom',     "pois",      NA),
#'                    nrow = s, ncol = s, byrow = TRUE)
#' f_dist <- array(f_dist_1, dim = c(s, s, degree + 1))
#' dim(f_dist)
#'
#' # ---------------------------------------------------------------------------
#' # Both p and f are drifting - Model 1.
#' # ---------------------------------------------------------------------------
#'
#' obj_fit_parametric <- fit_dsmm(sequence = sequence,
#'                                states = states,
#'                                degree = degree,
#'                                f_is_drifting = TRUE,
#'                                p_is_drifting = TRUE,
#'                                initial_dist = 'unif',
#'                                estimation = 'parametric',
#'                                f_dist = f_dist)
#' cat("The class of `obj_fit_parametric` is : (",
#'     paste0(class(obj_fit_parametric), collapse = ', '), ").\n")
#' # Estimated parameters.
#' f_params <- obj_fit_parametric$dist$f_drift_parameters
#' # The drifting sojourn time distribution parameters.
#' f_0 <- f_params[,,,1]
#' f_1.3 <- f_params[,,,2]
#' f_2.3 <- f_params[,,,3]
#' f_1 <- f_params[,,,4]
#'
#'
#' params <- paste0('q = ', round(f_params["c", "t", 1, ], 3),
#'                  ', beta = ', round(f_params["c", "t", 2, ], 3))
#' f_names <- c("f_0", paste0("f_", 1:(degree-1), "/", degree), "f_1")
#' all_names <- paste(f_names, ":", params)
#' cat("The drifting of the parameters for passing from \n",
#'     "`u` = 'c' to `v` = 't' under a discrete Weibull distribution is:",
#'     "\n", all_names[1], "\n", all_names[2],
#'     "\n", all_names[3], "\n", all_names[4])
#'
#' # ---------------------------------------------------------------------------
#' # f is not drifting, only p is drifting - Model 2.
#' # ---------------------------------------------------------------------------
#'
#' obj_fit_parametric_2 <- fit_dsmm(sequence = sequence,
#'                                states = states,
#'                                degree = degree,
#'                                f_is_drifting = FALSE,
#'                                p_is_drifting = TRUE,
#'                                initial_dist = 'unif',
#'                                estimation = 'parametric',
#'                                f_dist = f_dist_1)
#' cat("The class of `obj_fit_parametric_2` is : (",
#'     paste0(class(obj_fit_parametric_2), collapse = ', '), ").\n")
#' # Estimated parameters.
#' f_params_2 <- obj_fit_parametric_2$dist$f_notdrift_parameters
#'
#'
#' params_2 <- paste0('q = ', round(f_params_2["c", "t", 1], 3),
#'                  ', beta = ', round(f_params_2["c", "t", 2], 3))
#'
#' cat("Not-drifting parameters for passing from ",
#'     "`u` = 'c' to `v` = 't' \n under a discrete Weibull distribution is:\n",
#'     paste("f :", params_2))
#'
#'
#' # ===========================================================================
#' # Some methods for the `dsmm_fit_nonparametric` and
#' #  `dsmm_fit_parametric` objects.
#' # ===========================================================================
#' sim_seq_nonparametric <- simulate(obj_model_1, nsim = 10)
#' str(sim_seq_nonparametric)
#'
#' kernel_drift_parametric <- get_kernel(obj_fit_parametric, klim = 10)
#' str(kernel_drift_parametric)
fit_dsmm <- function(sequence, states, degree,
                     f_is_drifting, p_is_drifting,
                     initial_dist = "unif",
                     estimation = "nonparametric",
                     f_dist = NULL
                     # , numerical_est = FALSE
                     ) {
    # Sequence
    if (missing(sequence)) {
        stop("\nPlease provide the `sequence` character ",
             "vector in order to fit the model.")
    }
    # State Space
    if (missing(states)) {
        stop("\nPlease provide a character vector for the",
    "`states` of the State Space.")
    }
    # Degree
    if (missing(degree)) {
        stop("\nPlease provide the polynomial `degree` as an integer.")
    }
    # f_is_drifting
    if (missing(f_is_drifting)) {
        stop("\nPlease provide whether the sojourn time distribution",
             "f is drifting through the logical parameter `f_is_drifting`.")
    } else if (!is_logical(f_is_drifting)) {
        stop("\nThe logical parameter `f_is_drifting` ",
             "should be either TRUE or FALSE.")
    }
    # p_is_drifting
    if (missing(p_is_drifting)) {
        stop("\nPlease provide whether the transition matrix p is drifting",
             " through the logical parameter `p_is_drifting`.")
    } else if (!is_logical(p_is_drifting)) {
        stop("\nThe logical parameter `p_is_drifting` ",
             "should be either TRUE or FALSE.")
    }
    ## ## This will be implemented in a future version of the package. ## ##
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ## if (missing(numerical_est)) {
    ##     stop("\nPlease provide whether the numarical estimation should be",
    ##         "used to fit the model, through the logical parameter ",
    ##         "`numerical_est`",
    ##         ".\nThis should be used only when f or p is drifting.")
    ## }
    ## if (!is_logical(numerical_est)) {
    ##     stop("\nThe logical parameter `numerical_est` should be either",
    ##          " TRUE or FALSE.")
    ## } else if (numerical_est == TRUE) {
    ##     stop("\nThe numerical estimation under constraint, ",
    ##          "currently defined by the ",
    ##          "argument `numerical_est` == TRUE, ",
    ##          "has not been implemented yet.")
    ## }
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    # Check for validity of the parameters to be used to be used.
    if (valid_degree(degree)) {
        degree <- as.integer(degree)
    }
    stopifnot(valid_states(states),
              valid_initial_dist_char(initial_dist),
              valid_model(p_is_drifting, f_is_drifting))
    # States
    states <- sort(states)
    s <- length(states)
    stopifnot(valid_sequence(sequence, s))
    # Estimation
    stopifnot(valid_estimation(estimation = estimation,
                               fpar = f_dist, s = s,
                               f_is_drifting = f_is_drifting,
                               degree = degree, states = states))
    # Degree
    D <- degree + 1L
    # Sequence
    seq_encoding <- rle(sequence)
    # Model Size = length(sequence) - 1
    n <- (N <- length(seq <- seq_encoding$values)) - 1
    # k_max
    k_max <- max(X <- seq_encoding$lengths)
    # Check for sojourn times all equal to an integer.
    if (length(X_unique <- unique(X)) == 1) {
        stop("\nPlease input a sequence with Sojourn times",
             " X not all equal to `", X_unique, "`.")
    }
    # Initial Distribution.
    if (initial_dist == "unif") {
        alpha <- rep(1/s, s)
    } else if (initial_dist == "freq") {
        alpha <- as.vector(table(seq) / N)
    }
    # # if we have many sequences... another way is possible:
    # alpha <- sapply(sapply(states, function(state) which(state == seq[1])),
    #                 function(which_states) as.double(length(which_states)))
    names(alpha) <- states
    # Model
    model <- get_model(p_is_drifting, f_is_drifting)
    # Estimation
    # Get `vector_1_u`, `vector_1_uvl`.
    vector_1_u <- get_1_u(seq, N, n, states, s)
    vector_1_uvl <- get_1_uvl(seq, N, n, X, k_max, states, s)
    # Get `A_i` and the matrix of multiplications, `Aij`.
    Ai <- get_A_i(degree, model_size = n)
    Aij <- apply(Ai, c(1), function(a_i) a_i * t(Ai))
    dim(Aij) <- c(N, D, D)
    # Get Mij according to theory.
    Mij <- apply(Aij[-1,,], c(2,3), function(a_ij) a_ij %*% vector_1_u)
    # dimnames(Mij) <- list(
    #     as.list(states),
    #     as.list(names_i_d(degree, 'A')),
    #     as.list(names_i_d(degree, 'A'))
    # )
    # Mij <- aperm(Mij, c(2,3,1)) # This is how Mij appears in theory...
    det_Mij <- apply(Mij, c(1), det)
    min_Mij <- min(det_Mij)
    if (all_equal_numeric(min_Mij, 0)) {
        stop("The system MJ = P does not have any solution for the state(s):",
             paste0('\n"', states[which(det_Mij <= 0)], '"'))
    }
    # Get matrix `P_i`
    Pi <- apply(Ai[,-1], c(1), function(a_i) a_i %*% vector_1_uvl)
    dim(Pi) <- c(s, s, k_max, D)
    # dimnames(Pi) <- list(
    #     as.list(paste0('u = ', states)),
    #     as.list(paste0('v = ', states)),
    #     as.list(paste0('l = ', 1:k_max)),
    #     as.list(names_i_d(degree, 'A'))
    # )
    Ji <- sapply(1:s,
                 function(u) apply(Pi[u,,,], c(1,2),
                                   function(Pi_uvl) solve(Mij[u,,], Pi_uvl)))
    dim(Ji) <- c(D, s, k_max, s)
    # dimnames(Ji) <- list(
    #     as.list(names_i_d(degree)),
    #     as.list(paste0('v = ', states)),
    #     as.list(paste0('l = ', 1:k_max)),
    #     as.list(paste0('u = ', states))
    # )
    Ji <- aperm(Ji, c(4, 2, 3, 1)) # Just for clarification.
    # Normalizing Values.
    Jinormal <- apply(Ji, c(1,4), function(vl_values) {
        if (any(vl_values < 0)) {
            return(get_f(vl_values))
        } else {
            return(vl_values)
        }
    })
    # dimnames(Jinormal) <- list(
    #     as.list(paste0("v:", states, ", l = ", rep(1:k_max, each = s))),
    #     as.list(states),
    #     as.list(names_i_d(degree, 'q'))
    # )
    dim(Jinormal) <- c(s, k_max, s, D)
    Ji <- aperm(Jinormal, c(3,1,2,4))
    dimnames(Ji) <- list(
        as.list(states),
        as.list(states),
        as.list(paste0("l = ", 1:k_max)),
        as.list(names_i_d(degree, 'q'))
    )
    # Return with `dist`
    p_drift <- apply(Ji, c(1,2,4), sum)
    dimnames(p_drift) <- NULL
    dimnames(p_drift) <- list(lstates <- as.list(states),
                              lstates, as.list(names_i_d(degree, 'p')))
    p_drift_array <- array(apply(p_drift, c(3),
                                 function(matrix_uv)
                                     rep(matrix_uv, times = k_max)),
                           dim = c(s, s, k_max, D))
    f_drift <- Ji / p_drift_array
    f_drift[is.nan(f_drift)] <- 0
    # Give names to the returned object.
    dimnames(f_drift) <- NULL
    dimnames(f_drift) <- list(lstates, lstates,
                              paste0('l = ', 1:k_max),
                              as.list(names_i_d(degree, 'f')))
    dist <- list('p_drift' = p_drift, 'f_drift' = f_drift)
    # END MODEL 1
    if (model == "Model_2") {
        # p_drift is already computed from Model 1!
        # Compute f_notdrift.
        f_numerator <- apply(Ji, c(1,2,3), sum) # sums over degree + 1.
        f_denominator <- apply(p_drift, c(1,2), sum) # sums over degree + 1.
        f_notdrift <- f_numerator / rep(f_denominator, k_max)
        f_notdrift[is.nan(f_notdrift)] <- 0
        # Return `dist`.
        dist <- list("p_drift" = p_drift, "f_notdrift" = f_notdrift)
        # printf('fnotdrift sums', apply(f_notdrift, c(1,2), sum))
        # END MODEL 2.
    } else if (model == "Model_3") {
        # BEGIN MODEL 3
        p_notdrift <- apply(p_drift, c(1,2), sum) / D
        f_drift_model3 <- Ji / c(p_notdrift)
        f_drift_model3[is.nan(f_drift_model3)] <- 0
        apply(f_drift_model3, c(1,2,4), sum)
        # Return `dist`.
        dist <- list("p_notdrift" = p_notdrift, "f_drift" = f_drift)
        # printf('pnotdrift sum', rowSums(p_notdrift))
        # END MODEL 3.
    }
    # Parametric Estimation.
    if (estimation == 'parametric') {
        if (f_is_drifting) {
            f_drift <- dist$f_drift
            f_param <- sapply(1:D, function(d) {
                f_ijl <- f_drift[,,,d]
                dist_ij <- f_dist[,,d]
                sapply(1:s, function(j) {
                    f_il <- f_ijl[,j,]
                    dist_i <- dist_ij[,j]
                    sapply(1:s, function(i) {
                            if (i == j) {
                                return(c(NA, NA)) # Skip diagonal values.
                            }
                            f_l <- f_il[i,]
                            par_dist <- dist_i[i]
                            sum_fl <- sum(f_l)
                            if (is.na(par_dist)) {
                                if (sum_fl != 0) {
                                    stop("For the states u = ", i, ", v = ", j,
                                         " and ", names_i_d(degree, 'f')[d],
                                         " `f_dist` was defined to be NA, ",
                                         "however we do have occurences in ",
                                         "the sequence, corresponding to ",
                                         "these states. Therefore, a ",
                                         "distribution should be chosen.")
                                }
                                return(c(NA, NA))
                            } else if (sum_fl == 0) {
                                stop("For the states u = ", i, ", v = ", j,
                                     " and ", names_i_d(degree, 'f')[d],
                                     " we have no occurences in the sequence.",
                                     " Therefore, `f_dist` has to be NA.")
                            }
                            suppressWarnings( # because NA is being generated.
                                parametric_estimation(
                                    lprobs = f_l, dist = par_dist, kmax = k_max,
                                    i = i, j = j, d = d, degree = degree,
                                    states = states
                                )
                            )
                    })
                })
            })
            dim(f_param) <- c(2, s, s, D)
            f_param <- aperm(f_param, c(2, 3, 1, 4))
            dimnames(f_param) <- list(
                lstates, lstates,
                paste0(1:2),
                names_i_d(degree, 'fpars')
            )
        } else {
            # not drifting case
            f_notdrift <- dist$f_notdrift
            f_param <- sapply(1:s, function(j) {
                f_il <- f_notdrift[,j,]
                dist_i <- f_dist[,j]
                sapply(1:s, function(i) {
                    if (i == j) {
                        return(c(NA, NA)) # Skip diagonal values.
                    }
                    f_l <- f_il[i,]
                    dist <- dist_i[i]
                    sum_fl <- sum(f_l)
                    if (is.na(dist)) {
                        if (sum_fl != 0) {
                            stop("For the states u = ", i, ", v = ", j,
                                 " and ", names_i_d(degree, 'f')[0],
                                 " `f_dist` was defined to be NA, however",
                                 " we do have occurences in the sequence.",
                                 " Therefore, a distribution should be",
                                 " chosen.")
                        }
                        return(c(NA, NA))
                    } else if (sum_fl == 0) {
                        stop("For the states u = ", i, ", v = ", j,
                             " and ", names_i_d(degree, 'f')[0],
                             " we have no occurences in the sequence.",
                             " Therefore, `f_dist` has to be NA.")
                    }
                    parametric_estimation(
                        lprobs = f_l, dist = dist, kmax = k_max,
                        i = i, j = j, d = 0,
                        degree = degree, states = states
                    )
                })
            })
            dim(f_param) <- c(2, s, s)
            f_param <- aperm(f_param, c(2, 3, 1))
            dimnames(f_param) <- list(
                lstates, lstates,
                paste0(1:2)
            )
        }
        # Changing the returned value.
        dist[[2]] <- f_dist
        dist$f_param <- f_param
        if (model == "Model_1") {
            names(dist)[2:3] <- c('f_drift_parametric', 'f_drift_parameters')
        } else if (model == "Model_2") {
            names(dist)[2:3] <- c('f_notdrift_parametric', 'f_notdrift_parameters')
        } else if (model == "Model_3") {
            names(dist)[2:3] <- c('f_drift_parametric', 'f_drift_parameters')
        }
    }
    # Get the `dist` parameter.
    obj <- list(
        "dist" = dist,
        "seq" = seq,
        "soj_times" = X,
        "k_max" = k_max,
        "initial_dist" = alpha,
        "model_size" = n,
        "states" = states,
        "s" = s,
        "degree" = degree,
        # "numerical_est" = numerical_est,
        "f_is_drifting" = f_is_drifting,
        "p_is_drifting" = p_is_drifting,
        "Model" = model,
        "estimation" = estimation,
        "A_i" = get_A_i(degree, model_size = n), # used in `get_kernel()`.
        "Ji" = Ji # The semi-markov kernels used for the estimations.
    )
    if (estimation == 'nonparametric') {
        class(obj) <- c("dsmm_fit_nonparametric", "dsmm")
    } else {
        class(obj) <- c("dsmm_fit_parametric", "dsmm")
    }
    return(obj)
}



