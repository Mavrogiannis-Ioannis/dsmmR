# '''
#    1. This file concerns itself with the creation and definition of the
#    fitted drifting semi-Markov model on a random sequence, called `dsmm_fit`.
#    It is a child class of `dsmm`.
#    It inherits the methods:
#    `simulate`, `get_kernel`.
#    It has its own methods for:
#    `print`, `check_attributes`.
# '''


#' @title Least Square Estimation (LSE) of a Drifting Semi-Markov Chain
#' @aliases dsmm_fit
#' @description Least Square Estimation of a Drifting semi-Markov chain,
#'     given one sequence of states. This estimation is non-parametric and
#'     the estimation of three different models is available.
#'
#' @param sequence Character vector that represents a sequence of states in
#'     \eqn{E}. States must be characters with \eqn{length \geq 1}.
#' @param states Character vector that represents the state space
#'     \eqn{E} of choice. It has length equal to \eqn{s}.
#' @param degree Positive integer that represents the polynomial degree
#'     \eqn{d} for the Drifting Semi-Markov Model.
#' @param f_is_drifting Logical. Specifies if \eqn{f} is drifting or not.
#' @param p_is_drifting Logical. Specifies if \eqn{p} is drifting or not.
#' @param initial_dist Optional. Character that describes the method to
#'     estimate the initial distribution.
#'     \itemize{
#'         \item \code{"unif"} : the initial probability of each state is
#'         equal to \eqn{1/s}.
#'         \item \code{"freq"} : the initial distribution of each state is
#'         equal to the frequencies of each state in the sequence.
#'     }
#' @param alt_est Logical. Specifies if an alternative estimation
#'     should be used for Models 2 and 3. \emph{Currently not supported}.
#'     Default value is \code{FALSE}.
#'
#' @details This function estimates a Drifting Semi-Markov Model the
#'     non-parametric case. Three possible models are able to be estimated.
#'     More about the three Models in \link{dsmmR}. A normalization technique
#'     is used in order to correct the estimation errors from small sequences.
#'
#' \strong{Model 1}
#'
#' When the Transition Matrix of the embedded Markov chain \eqn{p} and
#' the Conditional Sojourn Time Distribution \eqn{f} are both drifting,
#' the Drifting Semi-Markov Kernel can \strong{estimated} as:
#' \deqn{\hat{q}_{\frac{t}{n}}(u,v,l) =
#' \sum_{i = 0}^{d}A_{i}\hat{q}_{\frac{i}{d}}(u,v,l),}
#' \eqn{\forall t \in \{0,\dots,n\}, \forall u,v\in E,
#' \forall l\in \{0,\dots, k_{max} \}},
#' where \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree
#' \eqn{d} (see \link{dsmmR}). The kernels
#' \eqn{\hat{q}_{\frac{i}{d}}(u,v,l), i = 0, \dots, d},
#' are estimated through LSE, and are obtained after solving the system:
#' \deqn{MJ = P,}
#' where
#' \itemize{
#' \item \eqn{M = (M_{ij})_{i,j \in E} =
#'     (\sum_{t=1}^{n}1_{u}(t)A_{i}A_{j})_{i,j \in E}};
#' \item \eqn{J = (J_i)_{i\in E}= (\hat{q}_{\frac{i}{d}}(u,v,l))_{i\in E}};
#' \item \eqn{P=(P_i)_{i\in E}=(\sum_{t=1}^{n}1_{uvl}(t)A_{i}(t))_{i\in E}},
#' }
#' where \eqn{1_{u}(t)} is the identifier function that has the value 1 if
#' at the instance \eqn{t} the previous state on the
#' instance \eqn{t-1} was \eqn{u}, and the value 0 otherwise.
#' \eqn{1_{uvl}(t)} is defined similarly and it has the value 1, if at
#' the instance \eqn{t} the current state is equal to \eqn{v},
#' the previous state on the
#' instance \eqn{t-1} is equal to \eqn{u} and a sojourn time of \eqn{l}
#' was spent on \eqn{u} before jumping to \eqn{v}. It has the value 0 in all
#' the other cases.
#'
#' In order to obtain the estimations of \eqn{\hat{p}_{\frac{i}{d}}(u,v)}
#' and \eqn{\hat{f}_{\frac{i}{d}}(u,v,l)}, so that the following is satisfied,
#' \deqn{\hat{q}_{\frac{i}{d}}(u,v,l) =
#'     \hat{p}_{\frac{i}{d}}(u,v)\hat{f}_{\frac{i}{d}}(u,v,l),}
#' \eqn{\forall u,v \in E, \forall l \in \{0,\dots, k_{max} \}}, where
#' \eqn{k_{max}} is the maximum realized sojourn time,
#' we use the following estimations:
#'
#' \deqn{\hat{p}_{\frac{i}{d}}(u,v) =
#'     \sum_{l = 0}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l),}
#' \deqn{\hat{f}_{\frac{i}{d}}(u,v,l) =
#'     \frac{\hat{q}_{\frac{i}{d}}(u,v,l)}{
#'          \sum_{l = 0}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l)}.}
#'
#'
#' \strong{Model 2}
#'
#' In this case, \eqn{p} is drifting and \eqn{f} is not drifting. Therefore,
#' the estimated Drifting Semi-Markov Kernel will be given by:
#' \deqn{\hat{q}_{\frac{t}{n}}^{(2)}(u,v,l) =
#' \sum_{i=0}^{d}\hat{q}_{\frac{i}{d}}^{(2)}(u,v,l),}
#' \eqn{\forall t \in \{0,\dots,n\}, \forall u,v\in E,
#' \forall l\in \{0,\dots, k_{max} \}},
#' where \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree
#' \eqn{d} (see \link{dsmmR}). Since \eqn{p} is drifting,
#' we define the estimation of \eqn{p} the same way as we did in Model 1,
#' and since \eqn{f} is \strong{not drifting}, we will symbolize it by
#' \eqn{f_{notdrift}}. Therefore,
#' \eqn{\forall u,v \in E, \forall l \in \{0,\dots, k_{max} \}},
#' we have the following estimations:
#'
#' \deqn{\hat{p}_{\frac{i}{d}}(u,v) =
#'     \sum_{l = 0}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l),}
#' \deqn{\hat{f}_{notdrift}(u,v,l) =
#'     \frac{\hat{q}_{\frac{i}{d}}(u,v,l)}{
#'        \sum_{i = 0}^{d}\sum_{l=0}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l)},}
#'
#' Thus, the \emph{estimated} kernels for Model 2,
#' \eqn{\hat{q}_{\frac{i}{d}}^{(2)}(u,v,l) =
#' \hat{p}_{\frac{i}{d}}(u,v)\hat{f}_{notdrift}(u,v,l)}, can be written with
#' regards to the \emph{estimated} kernels of Model 1,
#' \eqn{\hat{q}_{\frac{i}{d}}}, as in the following:
#'
#' \deqn{\hat{q}_{\frac{i}{d}}^{(2)}(u,v,l) = \frac{
#' \sum_{l=1}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l)
#' \sum_{i = 0}^{d}\hat{q}_{\frac{i}{d}}(u,v,l)}{
#' \sum_{i = 0}^{d}\sum_{l=1}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l)}.}
#'
#'
#' \strong{Model 3}
#'
#' In this case, \eqn{f} is drifting and \eqn{p} is not drifting. Therefore,
#' the estimated Drifting Semi-Markov Kernel will be given by:
#' \deqn{\hat{q}_{\frac{t}{n}}^{(3)}(u,v,l) =
#' \sum_{i=0}^{d}\hat{q}_{\frac{i}{d}}^{(3)}(u,v,l),}
#' \eqn{\forall t \in \{0,\dots,n\}, \forall u,v\in E,
#' \forall l\in \{0,\dots, k_{max} \}},
#' where \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree
#' \eqn{d} (see \link{dsmmR}). Since \eqn{f} is drifting,
#' we define the estimation of \eqn{f} the same way as we did in Model 1,
#' and since \eqn{p} is \strong{not drifting}, we will symbolize it by
#' \eqn{p_{notdrift}}. Therefore,
#' \eqn{\forall u,v \in E, \forall l \in \{0,\dots, k_{max} \}},
#' we have the following estimations:
#'
#' \deqn{\hat{p}_{notdrift}(u,v) =
#' \frac{\sum_{i=0}^{d}\sum_{l=0}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l)}{d+1}}
#' \deqn{\hat{f}_{\frac{i}{d}}(u,v,l) =
#'     \frac{\hat{q}_{\frac{i}{d}}(u,v,l)}{
#'          \sum_{l = 0}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l)}.}
#'
#' Thus, the \emph{estimated} kernels for Model 3,
#' \eqn{\hat{q}_{\frac{i}{d}}^{(3)}(u,v,l) =
#' \hat{p}_{notdrift}(u,v)\hat{f}_{\frac{i}{d}}(u,v,l)}, can be written with
#' regards to the \emph{estimated} kernels of Model 1,
#' \eqn{\hat{q}_{\frac{i}{d}}}, as in the following:
#'
#' \deqn{\hat{q}_{\frac{i}{d}}^{(3)}(u,v,l) = \frac{
#' \hat{q}_{\frac{i}{d}}(u,v,l)
#' \sum_{i=0}^{d}\sum_{l=0}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l)}
#' {(d+1)\sum_{l=1}^{k_{max}}\hat{q}_{\frac{i}{d}}(u,v,l)}.}
#'
#'
#'
#' @return Returns an object of S3 class \code{`dsmm_fit`, `dsmmm`}.
#' It has the following attributes:
#' \itemize{
#' \item \code{dist} : List that contains the \emph{p} and \emph{f}
#' estimated distributions from the sequence, given in the argument
#' \code{`sequence`};
#' \item \code{seq} : Character vector that contains the
#' \strong{state jumps of the original sequence}. It is this attribute of the
#' object that describes the length of the model \eqn{n}. Last state is also
#' included, for a total length of \eqn{n+1}, but it should not used;
#' \item \code{soj_times} : Numerical vector that contains the sojourn times
#' spent for each state in \code{`seq`} before the jump to the next state;
#' Last state is also
#' included, for a total length of \eqn{n+1}, but it should not used;
#' \item \code{k_max} : Numerical value that contains the maximum sojourn
#' time, meaning the maximum value in \code{`soj_times`};
#' \item \code{initial_dist} : Numerical vector that contains an estimation
#' for the initial distribution of the realized states in \code{`sequence`};
#' \item \code{model_size} : Integer value that contains the length of the
#' model This is equal to \code{length(seq) - 1}, for \code{`seq`} as defined
#' above. \item \code{states} : Character vector that contains the realized
#' states given in the argument \code{`sequence`};
#' \item \code{s} : Integer that contains the length of the state
#' space \eqn{E} given in the attribute \code{`states`}.
#' \item \code{degree} : Integer that contains the polynomial degree
#' \eqn{d} considered for the drifting of the model.
#' \item \code{f_is_drifting} : Logical. Passing down from the arguments.
#' \item \code{p_is_drifting} : Logical. Passing down from the arguments.
#' \item \code{alt_est} : Logical. Passing down from the arguments.
#' \item \code{Model} : Character vector. Possible values:
#' \code{["Model 1", "Model 2", "Model 3"]}, corresponding to whether
#' \eqn{p,f} are drifting or not.
#' \item \code{A_i} : Numerical Matrix. Used for the methods defined for the
#' object. Is not printed when viewing the object.
#' \item \code{Ji} : Numerical Array. Used for the methods defined for the
#' object. Is not printed when viewing the object.
#' }
#'
#' @seealso
#' For sequence simulation: \link{simulate.dsmm} and \link{create_sequence}.
#'
#' For more theory regarding Drifting Semi-Markov Models:
#' \link{dsmmR}, \link{get_kernel}, \link{parametric_dsmm},
#' \link{nonparametric_dsmm}
#'
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov
#' Models Toward Applications - Their Use in Reliability and DNA Analysis.
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#'
#'
#'
#' @export
#'
#' @examples
#' # Alternatively, we can use `data("lambda", package = "dsmmR")`
#' # to retrieve a sequence.
#' sequence <- create_sequence("DNA", len = 1000L)
#' states <- sort(unique(sequence))
#' degree <- 3
#'
#' # Fitting the sequence with the *first* model ===============================
#' obj_model_1 <- fit_dsmm(sequence = sequence,
#'                        states = states,
#'                        degree = degree,
#'                        f_is_drifting = TRUE,
#'                        p_is_drifting = TRUE,
#'                        initial_dist = "freq",
#'                        alt_est = FALSE # this is the default value.
#'                        )
#' cat(paste0(
#'     "We fitted a sequence with ", obj_model_1$Model, ",\n",
#'     "model size: n = ", obj_model_1$model_size, ",\n",
#'     "length of state space: s = ", obj_model_1$s, ",\n",
#'     "maximum sojourn time: k_max = ", obj_model_1$k_max, " and\n",
#'     "polynomial (drifting) Degree: d = ", obj_model_1$degree, "\n"
#' ))
#' # Get the drifting p and f arrays.
#' p_drift <- obj_model_1$dist$p_drift
#' f_drift <- obj_model_1$dist$f_drift
#' cat(paste0(
#'     "Dimension of p_drift: (s, s, d + 1) = (",
#'     paste(dim(p_drift), collapse = ", "), ");\n",
#'     "Dimension of f_drift: (s, s, k_max, d + 1) = (",
#'     paste(dim(f_drift), collapse = ", "), ")"
#' ))
#'
#' # Fitting the sequence with the *second* model ==============================
#' obj_model_2 <- fit_dsmm(sequence = sequence,
#'                        states = states,
#'                        degree = degree,
#'                        initial_dist = "unif",
#'                        f_is_drifting = FALSE,
#'                        p_is_drifting = TRUE)
#' cat("We fitted a sequence with", obj_model_2$Model, "\n")
#' # Get the drifting p and non-drifting f arrays.
#' p_drift_2 <- obj_model_2$dist$p_drift
#' f_notdrift <- obj_model_2$dist$f_notdrift
#' all.equal.numeric(p_drift, p_drift_2) # p is the same as in Model 1.
#' cat(paste0(
#'     "Dimension of f_notdrift: (s, s, k_max) = (",
#'      paste(dim(f_notdrift), collapse = ", "), ")"
#' ))
#'
#' # Fitting the sequence with the *third* model ===============================
#' obj_model_3 <- fit_dsmm(sequence = sequence,
#'                        states = states,
#'                        degree = degree,
#'                        f_is_drifting = TRUE,
#'                        p_is_drifting = FALSE)
#' cat("We fitted a sequence with", obj_model_3$Model, "\n")
#' # Get the drifting f and non-drifting p arrays.
#' p_notdrift <- obj_model_3$dist$p_notdrift
#' f_drift_3 <- obj_model_3$dist$f_drift
#' all.equal.numeric(f_drift, f_drift_3) # f is the same as in Model 1.
#' cat(paste0(
#'     "Dimension of f_notdrift: (s, s) = (",
#'     paste(dim(p_notdrift), collapse = ", "), ")"
#' ))
#'
#' # Some methods for the `dsmm_fit` objects.
#' k1 <- get_kernel(obj_model_1)
#' str(k1)
#'
#' sim_seq <- simulate(obj_model_1, nsim = 50)
#' str(sim_seq)
fit_dsmm <- function(sequence, states, degree,
                     f_is_drifting, p_is_drifting,
                     initial_dist = "unif",
                     alt_est = FALSE) {
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
    # alt_methods
    ## ## This will be implemented in a future version of the package. ## ##
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ## if (missing(alt_est)) {
    ##     stop("\nPlease provide whether the alternative estimation should be",
    ##         "used to fit the model, through the logical parameter `alt_est`",
    ##         ".\nThis should be used only when f or p is drifting.")
    ## }
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    if (!is_logical(alt_est)) {
        stop("\nThe logical parameter `alt_est` should be either",
             " TRUE or FALSE.")
    } else if (isTRUE(alt_est)) {
        stop("\nThe alternative estimation, currently defined by the ",
             "argument `alt_est` = TRUE, has not yet been implemented.")
    }
    # Check for validity of the paremeters to be used to be used.
    stopifnot(valid_states(states),
              valid_degree(degree),
              valid_initial_dist_char(initial_dist),
              valid_model(p_is_drifting, f_is_drifting, alt_est))
    # Setup
    # States
    states <- sort(states)
    s <- length(states)
    stopifnot(valid_sequence(sequence, s))
    # Degree
    degree <- as.integer(degree)
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
    # Initial Distribution. (Perhaps it's slow...)
    if (initial_dist == "unif") {
        alpha <- sapply(sapply(states, function(state)
            which(state == seq[1])),
            function(which_states)
                as.double(length(which_states)))
    } else if (initial_dist == "freq") {
        alpha <- table(seq)
    }
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
    Jinormal <- apply(Ji, c(1,4), get_f)
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
        # END MODEL 2.
        # printf('fnotdrift sums', apply(f_notdrift, c(1,2), sum))
    } else if (model == "Model_3") {
        # BEGIN MODEL 3
        p_notdrift <- apply(p_drift, c(1,2), sum) / D
        f_drift_model3 <- Ji / c(p_notdrift)
        f_drift_model3[is.nan(f_drift_model3)] <- 0
        apply(f_drift_model3, c(1,2,4), sum)
        # Return `dist`.
        dist <- list("p_notdrift" = p_notdrift, "f_drift" = f_drift)
        # END MODEL 3.
        # printf('pnotdrift sum', rowSums(p_notdrift))
    }
    ## END ESTIMATION 1.
    # } else {
    #     ### BEGIN ESTIMATION 2.
    #     if (Model_2) {
    #         # BEGIN MODEL 2.
    #         # END MODEL 2.
    #     }
    #     if (Model_3) {
    #         ## BEGIN MODEL 3
    #         # END MODEL 3
    #     }
    #     ### END ESTIMATION 2.
    # }
    ##### Get the `dist` parameter.
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
        "f_is_drifting" = f_is_drifting,
        "p_is_drifting" = p_is_drifting,
        "alt_est" = alt_est,
        "Model" = model,
        "A_i" = get_A_i(degree, model_size = n), # used in `get_kernel()`.
        "Ji" = Ji
    )
    class(obj) <- c("dsmm_fit", "dsmm")
    return(obj)
}
