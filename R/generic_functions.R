# '''
#    1. This file contains the generic functions that are used for the objects
#    defined in this package.
#    These objects are of the class `dsmm`, which acts like the parent class,
#    and then the 3 child classes:
#       `dsmm_fit`, `dsmm_nonparametric` and `dsmm_parametric`.
#    The `dsmm` class is only used if there is no need to classify a function
#    for the three child classes.
#    2. It is worth noting that, from the following functions, ONLY the
#    generics have documentation available to the user,
#    apart from the functions :
#       `get_kernel()`, `simulate.dsmm()`, `is.dsmm()`,
#       `is.dsmm_fit()`, `is.dsmm_nonparametric()` and `is.dsmm_parametric()`,
#    which are all available to the user.
# '''


# ______________________________________________________________________________
# Checking the validity of the attributes passed into the functions.
# ______________________________________________________________________________
check_attributes <- function(obj) UseMethod('check_attributes', obj)
check_attributes.dsmm_fit <- function(obj) {
    # Check for names of the object.
    # Check whether `f_is_drifting`, `p_is_drifting` and
    # `alt_est`, are correctly used.
    if (!is_logical(f_is_drifting <- obj$f_is_drifting)) {
        stop("\nThe logical parameter `f_is_drifting` ",
             "should be either TRUE or FALSE.")
    }
    if (!is_logical(p_is_drifting <- obj$p_is_drifting)) {
        stop("\nThe logical parameter `p_is_drifting` ",
    "should be either TRUE or FALSE.")
    }
    # alt_est
    if (!is_logical(alt_est <- obj$alt_est)) {
        stop("\nThe logical parameter `alt_est` should ",
             "be either TRUE or FALSE.")
    }
    # Check names of `dist`, in order for `obj$dist[[ i ]]`
    # to work in the next section of `stopifnot`.
    pname <- if (p_is_drifting) 'p_drift' else 'p_notdrift'
    fname <- if (f_is_drifting) 'f_drift' else 'f_notdrift'
    if (any((names_tmp <- names(obj$dist)) !=
            (names_real <- c(pname, fname)))) {
        stop('\n`', paste0(substitute(obj), '$dist`'),
             ' should have the named distributions in order, as in: \n',
             paste0(1:2, ". ", names_real, collapse = ', '),
             ', when it has :\n',
             paste0(1:2, ". ", names_tmp, collapse = ', '))
    }
    stopifnot(
        valid_seq(seq = (seq <- obj$seq)),
        valid_model_size(model_size = obj$model_size,
                         length_seq = (l <- length(seq))),
        valid_soj_times(soj_times = (soj_times <- obj$soj_times),
                        length_seq = l),
        valid_k_max(k_max = obj$k_max, soj_times = soj_times),
        valid_states(states = (states <- obj$states)),
        valid_length_states(s = obj$s, states = states),
        valid_degree(degree = obj$degree),
        valid_model(p_is_drifting = p_is_drifting,
                    f_is_drifting = f_is_drifting,
                    alt_est = alt_est)
    )
    TRUE
}
check_attributes.dsmm_nonparametric <- function(obj) {
    # Check for names of the object.
    if (!is_logical(f_is_drifting <- obj$f_is_drifting)) {
        stop("\nThe logical parameter `f_is_drifting` should be ",
             "either TRUE or FALSE.")
    }
    if (!is_logical(p_is_drifting <- obj$p_is_drifting)) {
        stop("\nThe logical parameter `p_is_drifting` should be ",
             "either TRUE or FALSE.")
    }
    # Check names of `dist`, in order for `obj$dist[[ i ]]`
    # to work in the next section of `stopifnot`.
    pname <- if (p_is_drifting) 'p_drift' else 'p_notdrift'
    fname <- if (f_is_drifting) 'f_drift' else 'f_notdrift'
    if (any((names_tmp <- names(obj$dist)) !=
             (names_real <- c(pname, fname)))) {
        stop('\n`', paste0(substitute(obj), '$dist`'),
             ' should have the named distributions in order, as in: \n',
             paste0(1:2, ". ", names_real, collapse = ', '),
             ', when it has :\n',
             paste0(1:2, ". ", names_tmp, collapse = ', '))
    }
    # Check whether the object attributes are correctly used.
    stopifnot(
        is_integer(obj$model_size),
        is_integer(k_max <- obj$k_max),
        valid_states(states <- obj$states),
        valid_length_states(s <- obj$s, states),
        valid_degree(degree <- obj$degree),
        valid_model(p_is_drifting, f_is_drifting),
        valid_p_dist(p_dist = obj$dist[[1]],
                     states = states,
                     s = s, degree = degree,
                     p_is_drifting = p_is_drifting),
        valid_fdist_nonparametric(f_dist = obj$dist[[2]],
                                  states = states, s = s,
                                  degree = degree,
                                  f_is_drifting = f_is_drifting,
                                  k_max = k_max)
    )
    TRUE
}

check_attributes.dsmm_parametric <- function(obj) {
    if (!is_logical(f_is_drifting <- obj$f_is_drifting)) {
        stop("\nThe logical parameter `f_is_drifting` should be either ",
             "TRUE or FALSE.")
    }
    if (!is_logical(p_is_drifting <- obj$p_is_drifting)) {
        stop("\nThe logical parameter `p_is_drifting` should be either ",
             "TRUE or FALSE.")
    }
    # Check names of dist, in order for `obj$dist[[ i ]]` to work in
    # the next section of `stopifnot`.
    pname <- if (p_is_drifting) 'p_drift' else 'p_notdrift'
    fname <- if (f_is_drifting) 'f_drift_parametric' else 'f_notdrift'
    fparname <-
        if (f_is_drifting) 'f_drift_parameters' else 'f_notdrift_parameters'
    if (any((names_tmp <- names(obj$dist)) !=
            (names_real <- c(pname, fname, fparname)))) {
        stop('\n`', paste0(substitute(obj), '$dist`'),
             ' should have the named distributions in order, as in: \n',
             paste0(1:2, ". ", names_real, collapse = ', '),
             ', when it has :\n',
             paste0(1:2, ". ", names_tmp, collapse = ', '))
    }
    # Check whether `alt_est`, `f_is_drifting` and `p_is_drifting`
    # are correctly used.
    stopifnot(is_integer(obj = obj$model_size),
              valid_model(p_is_drifting, f_is_drifting),
              valid_states(states = (states <- obj$states)),
              valid_length_states(s = (s <- obj$s), states),
              valid_initial_dist(initial_dist = obj$initial_dist,
                                         s = s),
              valid_degree(degree = (degree <- obj$degree)),
              valid_p_dist(p_dist = obj$dist[[1]],
                           states = states,
                           s = s, degree = degree,
                           p_is_drifting = p_is_drifting),
              valid_fdist_parametric(fdist = obj$dist[[2]],
                                     params = obj$dist[[3]],
                                     s = s, degree = degree,
                                     f_is_drifting = f_is_drifting)
    )
    TRUE
}



# ______________________________________________________________________________
# Checking the class of an object and if it satisfies the necessary conditions.
# ______________________________________________________________________________

#' @title Function to check if an object is of class \code{dsmm}
#'
#' @description It checks for inheritance of the class \code{dsmm} and
#' also for the validity of the attributes specified. This class acts like
#' a parent class for the classes \code{dsmm_fit, dsmm_parametric,
#' dsmm_nonparametric}.
#'
#' @param obj Arbitrary \code{R} object.

#' @seealso \link{is.dsmm_fit}, \link{is.dsmm_parametric},
#' \link{is.dsmm_nonparametric}
#'
#' @return TRUE or FALSE.
#'
#' @export
is.dsmm <- function(obj) {
    # Check for missing object, `obj`.
    if (missing(obj)) {
        stop("\nPlease input the `obj` parameter of class `dsmm`. ",
             "This can be done through the functions\n `fit_dsmm`,",
             "`dsmm_parametric`")
    }
    # Check for class of object.
    if (!inherits(obj, 'dsmm')) {
        cat("\n`obj` is not of class `dsmm_fit`.")
        return()
    }
    check_attributes(obj)
}

#' @title Function to check if an object is of class \code{dsmm_fit}
#'
#' @description It checks for inheritance of the class \code{dsmm_fit} and
#' also for the validity of the attributes specified.
#'
#' @param obj Arbitrary \code{R} object.
#'
#' @seealso \link{is.dsmm}, \link{is.dsmm_parametric},
#' \link{is.dsmm_nonparametric}
#'
#' @return TRUE or FALSE.
#'
#' @export
is.dsmm_fit <- function(obj) {
    # Check for missing object, `obj`.
    if (missing(obj)) {
        stop("\nPlease input the `obj` parameter of class `dsmm_fit`.\n",
             "This can be done through the function `fit_dsmm`.")
    }
    # Check for class of object.
    if (!inherits(obj, 'dsmm_fit')) {
        cat("\n`obj` is not of class `dsmm_fit`.")
        return()
    }
    check_attributes(obj)
}

#' @title Function to check if an object is of class
#' \code{dsmm_nonparametric}
#'
#' @description It checks for inheritance of the class
#' \code{dsmm_nonparametric} and also for the validity of
#' the attributes specified. This class inherits methods from the parent class
#' \code{dsmm}.
#'
#' @seealso \link{is.dsmm}, \link{is.dsmm_fit}, \link{is.dsmm_parametric}
#'
#' @param obj Arbitrary \code{R} object.
#'
#' @return TRUE or FALSE.
#'
#' @export
is.dsmm_nonparametric <- function(obj) {
    # Check for missing object, `obj`.
    if (missing(obj)) {
        stop("\nPlease input the `obj` parameter of class `dsmm`,",
             "`fit_dsmm`, `dsmm_parametric`, or `dsmm_nonparametric`.")
    }
    # Check for class of object.
    if (!inherits(obj, c('dsmm_nonparametric', 'dsmm'))) {
        stop("\n`obj` needs to be of class `dsmm_nonparametric`.")
    }
    check_attributes(obj)
}

#' @title Function to check if an object is of class \code{dsmm_parametric}
#'
#' @description It checks for inheritance of the class \code{dsmm_parametric}
#' and also for the validity of the attributes specified.
#'
#' @param obj Arbitrary \code{R} object.
#'
#' @seealso \link{is.dsmm}, \link{is.dsmm_fit},
#' \link{is.dsmm_nonparametric}
#'
#' @return TRUE or FALSE.
#'
#' @export
is.dsmm_parametric <- function(obj) {
    if (missing(obj)) {
        stop("\nPlease input the `obj` parameter of class `dsmm_parametric`.",
             " This can be done through the function `dsmm_parametric`.")
    }
    # Check for class of object.
    if (!inherits(obj, c('dsmm_parametric', 'dsmm'))) {
        stop("\n`obj` needs to be of class `dsmm_parametric`. ",
             "This can be done through the function `dsmm_parametric`.")
    }
    check_attributes(obj)
}



# ______________________________________________________________________________
# Get the kernel q_(t/n) (u,v,l) that is necessary for the `simulate` function.
# ______________________________________________________________________________

#' @title Obtain the Drifting Semi-Markov kernel
#' @description
#' This is a generic method that computes and returns the Drifting
#' Semi-Markov kernel as a numerical array of dimensions (s, s, k_max, n + 1).
#'
#' @details
#' The Drifting Semi-Markov kernel \eqn{q_{\frac{t}{n}}(u,v,l)} describes
#' the probability to jump to the current state \eqn{v} when the previous
#' state \eqn{u} had a sojourn time of \eqn{l}, for every point \eqn{t} of a
#' model with length \eqn{n}.
#' Specifically, it is given as the sum of a linear combination:
#' \deqn{q_{\frac{t}{n}}(u,v,l)=
#'      \sum_{i = 0}^{d}A_{i}(t)q_{\frac{i}{d}}(u,v,l),}
#' where \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree
#' \eqn{d} that satisfy certain conditions (\link{dsmmR}) and
#' \eqn{q_{\frac{i}{d}}(u,v,l), i = 0, \dots, d}
#' are \eqn{d + 1} kernels that describe a multiplication between
#' a number of Markov Transition Matrices \eqn{p} and
#' a number of Sojourn Time Distributions \eqn{f}.
#' Three possible model specifications are described below.
#'
#' \strong{\emph{Model 1}}
#'
#' In this case, both \eqn{p} and \eqn{f} are "drifting" between \eqn{d + 1}
#' fixed points of the model, hence the "drifting" in
#' Drifting Semi-Markov Models. Therefore, we denote specifically this kernel,
#' resulting from \eqn{d + 1} \eqn{p} and \eqn{f} with \eqn{q_{\frac{t}{n}}^{(1)}} and
#' the corresponding Semi-Markov kernels with \eqn{q_{\frac{i}{d}}^{(1)}}.
#' \deqn{q_{\frac{i}{d}}^{(1)}(u,v,l)=
#'      {p_{\frac{i}{d}}(u,v)}{f_{\frac{i}{d}}(u,v,l)},}
#' where for \eqn{i = 0, \dots, d} we have \eqn{d + 1} Markov
#' Transition Matrices \eqn{p_{\frac{i}{d}}(u,v,l)},
#' and \eqn{d + 1} Sojourn Time Distributions
#' \eqn{f_{\frac{i}{d}}(u,v,l)}, where \eqn{d} is the polynomial degree.
#'
#' Thus, the Drifting Semi-Markov kernel will be equal to:
#'
#' \deqn{q_{\frac{t}{n}}^{(1)}(u,v,l) =
#'   \sum_{i = 0}^{d} A_i(t)q_{\frac{i}{d}}^{(1)}(u,v,l) =
#'   \sum_{i = 0}^{d} A_i(t)p_{\frac{i}{d}}(u,v)f_{\frac{i}{d}}(u,v,l)
#'    }
#'
#'
#' \strong{\emph{Model 2}}
#'
#' In this case, \eqn{p} is drifting and \eqn{f} \strong{is not drifting}.
#' Therefore, we denote specifically this kernel,
#' resulting from \eqn{d + 1} \eqn{p} and \eqn{f} with \eqn{q_{\frac{t}{n}}^{(2)}} and
#' the corresponding Semi-Markov kernels with \eqn{q_{\frac{i}{d}}^{(2)}}.
#' \deqn{q_{\frac{i}{d}}^{(2)}(u,v,l)={p_{\frac{i}{d}}(u,v)}{f(u,v,l)},}
#' where the \eqn{f(u,v,l)} is constant, not drifting alongside the point
#' \eqn{t} of the model with length \eqn{n}. Thus, the Drifting Semi-Markov kernel
#'  will be equal to:
#'
#' \deqn{q_{\frac{t}{n}}^{(2)}(u,v,l) =
#'   \sum_{i = 0}^{d} A_i(t)q_{\frac{i}{d}}^{(2)}(u,v,l) =
#'   \sum_{i = 0}^{d} A_i(t)p_{\frac{i}{d}}(u,v)f_{\frac{i}{d}}(u,v,l)
#'    }
#'
#'
#' \strong{\emph{Model 3}}
#'
#' In this case, \eqn{f} is drifting and \eqn{p} \strong{is not drifting}.
#' Therefore, the Markov kernel is now described as:
#' \deqn{q_{\frac{i}{d}}^{(3)}(u,v,l)={p(u,v)}{f_{\frac{i}{d}}(u,v,l)},}
#' where \eqn{p(u,v)} is constant, not drifting alongside the point
#' \eqn{t} of the model with length \eqn{n}. Thus, the Drifting Semi-Markov kernel
#'  will be equal to:
#'
#' \deqn{q_{\frac{t}{n}}^{(3)}(u,v,l) =
#'   \sum_{i = 0}^{d} A_i(t)q_{\frac{i}{d}}^{(3)}(u,v,l) =
#'   \sum_{i = 0}^{d} A_i(t)p_{\frac{i}{d}}(u,v)f_{\frac{i}{d}}(u,v,l)
#'    }
#'
#'
#' @param obj An object of S3 class \code{dsmm, dsmm_fit,
#' dsmm_nonparametric} or \code{dsmm_parametric}.
#' @param t Optional, but recommended. Positive integer specifying a
#' single point \eqn{t} in the length of the model.
#' @param u Optional. Character specifying the previous state \eqn{u}.
#' It can also be given as a positive integer, describing the order it appears
#' in the state space.
#' @param v Optional. Character specifying the next state \eqn{v}.
#' It can also be given as a positive integer, describing the order it appears
#' in the state space.
#' @param l Optional. Positive integer specifying the sojourn time \eqn{l}
#' that is spent in the previous state \eqn{u}.
#' @param klim Optional. Positive integer.
#' Used only for the S3 class \code{dsmm_parametric}.
#' Specifies the time horizon used to approximate the \eqn{d + 1}
#' (or 1, if \eqn{f} is \emph{not drifting}) sojourn time distributions.
#' Default value is 80L. A larger value will result in a considerably larger
#' kernel (which has dimensions of \eqn{(s, s, klim, n + 1)}), that will
#' require more memory and will slow down considerably the
#' \code{simulate.dsmm()} method.
#' (\link{dsmm_parametric}, \link{simulate.dsmm})
#'
#' @return An array with dimensions of \eqn{(n+1, s, s, k_{max})}, giving the
#' value of the Drifting Semi-Markov kernel \eqn{q_{\frac{t}{n}}(u,v,l)} for
#' the corresponding \eqn{(t,u,v,l)}. If any of \eqn{t,u,v,} or \eqn{l} were
#' specified, their dimension in the array becomes 1.
#'
#' @export
#'
#' @seealso
#' This kernel can be a result either from the estimation using: \link{fit_dsmm},
#' or through the functions defining a Drifting Semi-Markov model specification:
#' \link{parametric_dsmm}, \link{nonparametric_dsmm}.
#'
#' For sequence simulation through this kernel: \link{simulate.dsmm}.
#'
#' For the theoretical background of Drifting Semi-Markov models: \link{dsmmR}.
#'
#' @examples
#' # Setup.
#' states <- c("Rouen", "Bucharest", "Samos", "Aigio", "Marseille")
#' seq <- create_sequence(states, probs = c(0.3, 0.1, 0.1, 0.3, 0.2))
#' obj_model_2 <- fit_dsmm(
#'     sequence = seq,
#'     states = states,
#'     degree = 3,
#'     f_is_drifting = FALSE,
#'     p_is_drifting = TRUE
#' )
#'
#' # Get the kernel.
#' kernel_model_2 <- get_kernel(obj_model_2)
#' cat(paste0(
#'     "If no further arguments are made, kernel has dimensions for all ",
#'     "u, v, l, t:\n",
#'     "(s, s, k_max, n + 1) = (",
#'     paste(dim(kernel_model_2), collapse = ", "), ")"
#' ))
#'
#' # Specifying `t`.
#' kernel_model_2_t <- get_kernel(obj_model_2, t = 100)
#' # kernel_model_2[,,,t=100]
#' cat(paste0(
#'     "If we specify t, the kernel has dimensions for all the remaining ",
#'     "u, v, l:\n(s, s, k_max) = (",
#'     paste(dim(kernel_model_2_t), collapse = ", "), ")"
#' ))
#'
#' # Specifying `t` and `u`.
#' kernel_model_2_tu <- get_kernel(obj_model_2, t = 2, u = "Aigio")
#' # kernel_model_2["Aigio",,,t=100]
#' cat(paste0(
#'     "If we specify t and u, the kernel has dimensions for all the ",
#'     " remaining v, l:\n(s, k_max) = (",
#'     paste(dim(kernel_model_2_tu), collapse = ", "), ")"
#' ))
#'
#'
#' # Specifying `t`, `u` and `v`.
#' kernel_model_2_tuv <- get_kernel(obj_model_2, t = 3,
#'                                  u = "Rouen", v = "Bucharest")
#' # kernel_model_2["Aigio","Bucharest",,t=100]
#' cat(paste0(
#'     "If we specify t, u and v, the kernel has dimensions for all l:\n",
#'     "(k_max) = (", paste(length(kernel_model_2_tuv), collapse = ", "), ")"
#' ))
#'
#' # It is possible to ask for any valid combination of `u`, `v`, `l` and `t`.
get_kernel <- function(obj, t, u, v, l, klim = 80L) {
    # Check for missing values & the validity of the attributes given.
    if (missing(obj)) {
        stop("\nPlease provide the object `obj`.")
    }
    stopifnot(check_attributes(obj))
    n <- obj$model_size
    N <- n + 1
    k_max <- obj$k_max
    states <- obj$states
    s <- length(states)
    if (!missing(t)) {
        stopifnot(is_integer(t))
        if (t > N) {
            stop("\n`t` cannot be larger than n + 1, where n = ", n)
        }
    }
    if (!missing(u)) {
        if (is.character(u)) {
            stopifnot(valid_state(u, states))
        } else if (is_integer(u) && u > s) {
            stop("\nThe previous state `u` is specified as the numbered ",
                 "state `", u,
                 "` in the state space of total length s = ",
                 s, ".")
        }
        u <- states[which(states == u)]
    }
    if (!missing(v)) {
        if (is.character(v)) {
            stopifnot(valid_state(v, states))
        } else if (is_integer(v) && v > s) {
            stop("\nThe previous state `v` is specified as the numbered",
                 " state `", v,
                 "` in the state space of total length s = ", s, ".")
        }
        v <- states[which(states == v)]
    }
    if (!missing(l)){
        stopifnot(is_integer(l))
        if (l > k_max) {
            stop("\nThe maximum sojourn time specified, `", l,
                 "` cannot be more than ", k_max, ".")
        }
    }
    if (!is_integer(klim)) {
        stop("\nAttribute `klim` should be a positive integer.")
    } else if (klim > 200) {
        warning("\nAttribute `klim` = ", klim, " will result in",
                " considerable memory requirements for a larger model size.")
    }
    # Possible defined parameters are passed down to `UseMethod`.
    UseMethod(generic = 'get_kernel', object = obj)
}

#' @export
get_kernel.dsmm <- function(obj, t, u, v, l, ...) {
    # =========================================================================.
    # '''
    #    This function is valid for both `dsmm_fit` AND `dsmm_nonparametric`.
    # '''
    # =========================================================================.
    # Get the correct distributions `dist`.
    D <- obj$degree + 1L
    dist <- obj$dist
    n <- obj$model_size
    k_max <- obj$k_max
    states <- obj$states
    s <- length(states)
    Ai <- obj$A_i # this should be a matrix with dim(A_i) == c(d+1, n+1)
    Model <- obj$Model
    pdist <- dist[[1]]
    fdist <- dist[[2]]
    # Get J_i with regards to Model.
    if (Model == 'Model_1') {
        Ji <- fdist * c(apply(pdist, c(3),
                              function(M_uv) rep(M_uv, times = k_max)))
    } else if (Model == "Model_2") {
        # dim(fdist) == (s, s, k_max) --> (s, s, k_max, d + 1).
        f_vector <- rep(fdist, D)
        # dim(pdist) == (s, s, d+1) --> (s, s, k_max, d + 1)
        p_vector <- apply(pdist, MARGIN = c(3),
                          FUN = function(x) rep(x, k_max))
        Ji <- f_vector * p_vector
    } else if (Model == "Model_3") {
        # dim(pdist) == (s, s) --> (s, s, k_max, d+1)
        p_vector <- rep(pdist, k_max * D)
        Ji <- fdist * p_vector
    }
    # Get kernel q_(t/n) (u,v,l).
    dim(Ji) <- c(s*s*k_max, D)
    # dimnames(Ji) <-list(
    #   as.list(paste0(rep(E, s), sort(rep(E, s)),
    #    sort(rep(1:k_max, s * s)))),
    #   as.list(names_i_d(degree, "q"))
    # )
    kernel <- get_valid_kernel(Ji, Ai, s, n, k_max, states)
    kernel[u, v, l, t]
}

#' @export
get_kernel.dsmm_parametric <- function(obj, t, u, v, l, klim = 80L) {
    # '''
    #    This function returns the kernel q_(t/n) with variables (t, u, v, l)!
    # '''
    # Does not need missing() case, since there is only one argument
    # and it needs to be passed down from the generic function.
    stopifnot(check_attributes(obj))
    # Get the correct distributions `dist`.
    pdist <- obj$dist[[1]] # [u,v]
    fdist <- obj$dist[[2]] # [u,v,]
    fpar <- obj$dist[[3]]  # [u,v,,]
    #  A large klim will create an abnormaly large
    #  kernel and also will slow down the computation for the method
    # `simulate.dsmm()` considerably.
    f_vector <- get_fdist_parametric(fdist, fpar, klim)
    f_vector[which(f_vector < 1e-10)] <- 0
    p_vector <- rep(pdist, each = klim)
    Ji <- f_vector * p_vector
    degree <- obj$degree
    D <- degree + 1
    dim(Ji) <- c(klim,
                 s <- obj$s,
                 s,
                 D)
    Ji <- aperm(Ji, c(2, 3, 1, 4))
    dim(Ji) <- c(s * s * klim, D)
    Ai <- get_A_i(degree, n <- obj$model_size)
    kernel <- Ji %*% Ai
    dim(kernel) <- c(s, s, klim, n + 1)
    dimnames(kernel) <- list(
        as.list(states <- obj$states),
        as.list(states),
        as.list(paste0('l = ', 1:klim)),
        as.list(paste0('t = ', 0:n))
    )
    if (!is_prob(kernel)) stop("\nThe final kernel is not a probability.")
    kernel[u, v, l, t]
}


# ______________________________________________________________________________
# Printing methods for each of the child classes.
# ______________________________________________________________________________
#' @export
print.dsmm_fit <- function(x, ...) {
    # '''
    #    This function was made for the purpose of NOT printing certain
    #    parameters. For example, the polynomials `A_i` are too long
    #    and would obscure the printing of the object.
    # '''
    check_attributes(x)
    nm <- names(x)
    for (i in seq_along(nm)) {
        if (!nm[i] %in% c('A_i', 'dist', 'Ji')) {
            cat(paste0("\n$", nm[i], "\n"))
            if (nm[i] %in% c('k_max', 'model_size', 's',
                             'degree', 'states',
                             'f_is_drifting', 'p_is_drifting',
                             'alt_est', 'Model')) {
                cat(x[[i]], "\n")
            } else {
                print(x[[i]])
            }
        } else if (nm[i] == 'dist') {
            nd <- names(x$dist)
            for (j in seq_along(nd)) {
                cat(paste0("\n\n\n$dist$", nd[j], "\n\n"))
                print(x$dist[[j]])
            }
        }
    }
    cat("\n\nClass:", paste0(class(x), collapse = ", "))
}

#' @export
print.dsmm_nonparametric <- function(x, ...) {
    # '''
    #    This function was made for the purpose of NOT printing
    #    certain parameters. For example, the polynomials `A_i`
    #    are too long and would obscure the printing of the object.
    # '''
    check_attributes(x)
    nm <- names(x)
    for (i in seq_along(nm)) {
        if (!nm[i] %in% c('A_i', 'dist', 'Ji')) {
            cat(paste0("\n$", nm[i], "\n"))
            if (nm[i] %in% c('k_max', 'model_size', 's', 'degree', 'states',
                             'f_is_drifting', 'p_is_drifting', 'Model')) {
                cat(x[[i]], "\n")
            } else {
                print(x[[i]])
            }
        } else if (nm[i] == 'dist') {
            nd <- names(x$dist)
            for (j in seq_along(nd)) {
                cat(paste0("\n\n\n$dist$", nd[j], "\n\n"))
                print(x$dist[[j]])
            }
        }
    }
    cat("\n\nClass:", paste0(class(x), collapse = ", "))
}

#' @export
print.dsmm_parametric <- function(x, ...) {
    # '''
    #    This function was made for the purpose of NOT printing certain
    #    parameters. For example, the polynomials `A_i` are too long
    #    and would obscure the printing of the object.
    # '''
    check_attributes(x)
    nm <- names(x)
    for (i in seq_along(nm)) {
        if (!nm[i] %in% c('A_i', 'dist')) {
            cat(paste0("\n$", nm[i], "\n"))
            if (nm[i] %in% c('model_size', 's', 'degree', 'states',
                             'f_is_drifting', 'p_is_drifting', 'Model')) {
                cat(x[[i]], "\n")
            } else {
                print(x[[i]])
            }
        } else if (nm[i] == 'dist') {
            nd <- names(x$dist)
            for (j in seq_along(nd)) {
                cat(paste0("\n\n\n$dist$", nd[j], "\n\n"))
                print(x$dist[[j]])
            }
        } else if (nm[i] == 'f_distribution_parameters') {
            print(x$f_distribution_parameters)
        }
    }
    cat("\n\nClass:", paste0(class(x), collapse = ", "))
}

# ______________________________________________________________________________
# Simulate a sequence from any `dsmm` object.
# ______________________________________________________________________________
#' @title Simulate a sequence under a Drifting Semi-Markov kernel.
#' @description Generic function that simulates a number of states \code{nsim}
#' under the rule of a Drifting Semi-Markov kernel, which is retrieved from the
#' object \code{obj}, which in turn inherits from the S3 class \code{dsmm}.
#'
#' @param object An object of S3 class \code{dsmm,dsmm_fit}
#' \code{dsmm_nonparametric} or \code{dsmm_parametric}.
#' @param nsim Optional. A positive integer specifying the number of simulations
#' that will make up the sequence, with the maximum value being the model size
#' that is specified in \code{obj}. The model size is also the default value.
#' @param seq_length Optional. A positive integer that will ensure the simulated
#' sequence will not have a \emph{total length} greater than \code{seq_length}
#' (however, it is possible for the total length to be \emph{less} than
#' \code{seq_length}).
#' @param seed Optional. An integer specifying the initialization of the random
#' number generator.
#' @param klim Optional. Positive integer. Passed down to \code{get_kernel}
#' for the parametric object, with class \code{dsmm_parametric}.
#' @param ... Attributes passed down from the \code{simulate} method.
#' Currently not used.
#'
#' @seealso
#' About random number generation in R \code{\link[base:RNG]{RNG}}.
#'
#' Fitting a model through a sequence from this function: \link{fit_dsmm}.
#'
#' For the theoretical background of Drifting Semi-Markov Models: \link{dsmmR}.
#'
#' @return A character vector based on \code{nsim} simulations, with a
#' maximum length of \code{seq_length}.
#' @export
#' @examples
#' # Setup.
#' seq <- create_sequence("DNA", len = 1000L)
#' states <- sort(unique(seq))
#' d <- 1
#' obj_model_3 <- fit_dsmm(sequence = seq,
#'                         states = states,
#'                         degree = d,
#'                         f_is_drifting = TRUE,
#'                         p_is_drifting = FALSE)
#'
#' # Using the method `simulate.dsmm()`.
#' simulated_seq <- simulate(obj_model_3, seed = 1)
#' short_sim <- simulate(obj = obj_model_3, nsim = 50, seed = 1)
#' cut_sim <- simulate(obj = obj_model_3, seq_length = 50, seed = 1)
#' str(simulated_seq)
#' str(short_sim)
#' str(cut_sim)
simulate.dsmm <- function(object, nsim = NULL, seed = NULL,
                          seq_length = NULL, klim = 80L, ...) {
    # Parameters Setup.
    if (missing(object)) {
        stop("\nPlease provide an objectect of class `dsmm`.")
    } else if (!inherits(object, c('dsmm', 'dsmm_fit', 'dsmm_nonparemetric',
                                'dsmm_parametric'))) {
        stop("\nPlease provide an objectect of class `dsmm` to use for the",
             " function `simulate`.",
             "\nThe objectect can be created through the functions ",
             "`parametric`, `nonparametric` and `fitdsmm`.")
    }
    if (is_integer(nsim) && is_integer(seq_length)) {
        stop("\nPlease specify only one of `nsim` or `seq_length` for ",
             "the simulation.")
    }
    if (is.null(nsim)) {
        nsim <- object$model_size
    } else if (!is_integer(nsim)) {
        stop("\nThe number of simulations `nsim`",
             "needs to be a positive integer.")
    }
    if (!is.null(seq_length) && !is_integer(seq_length)) {
        stop("\nThe final length of the sequence `seq_length`",
             "needs to be a positive integer.")
    } else if (is_integer(seq_length)) {
        nsim <- seq_length
    }
    if (!is.null(seed)) {
        set.seed(seed)
    } # Otherwise, set.seed is automatic through user's `Sys.time()`.
    if (!is_integer(klim)) {
        stop("\nAttribute `klim` should be a positive integer.")
    } else if (klim > 100) {
        warning("\nAttribute `klim` = ", klim, " will result in",
                " considerable more computation time\n for the",
                " simulation of a larger sequence.")
    }
    stopifnot(check_attributes(object))
    # In order to get the first letter of the sequence.
    alpha <- object$initial_dist
    states <- object$states
    initial_state <- sample(states, prob = alpha, size = 1)
    s <- length(states)
    n <- object$model_size
    if (nsim > n) {
        stop("\nThe number of simulations `nsim` = ", nsim,
             "cannot be larger than the model size, n = ", n)
    }
    kernel <- get_kernel(object, klim = klim)
    k_max <- dim(kernel)[3] # We get `k_max` even for the parametric case.
    # Get the rest of the sequence.
    vl_names <- paste(states, sort(rep(1:k_max, s)))
    vl_vector <- c()
    u <- initial_state
    for (time in 1:nsim) { # max(nsim) == l!
        u <- states[which(states == u)]
        vl_vector <- c(vl_vector,
                       vl <- get_vl(t = time, u = u,
                                    vl_names, kernel, states))
        u <- vl[1]
    }
    l_vl <- length(vl_vector)
    seq <- c(initial_state, vl_vector[seq(1, l_vl, 2)])
    seq <- seq[-length(seq)]
    #remove last state... It does not have a corresponding sojourn time.
    X <- as.numeric(c(vl_vector[seq(2, l_vl, by = 2)]))
    new_seq <- unlist(sapply(seq_along(seq), function(i) rep(seq[i], X[i])))
    if (is.null(seq_length)) {
        return(new_seq)
    } else {
        return(new_seq[1:seq_length])
    }
}

