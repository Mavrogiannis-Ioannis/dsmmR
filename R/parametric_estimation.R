# '''
#    This file concerns itself with the helper functions for the parametric
#    estimation given a sequence. These are all used in the function `fit_dsmm`,
#    and the parametric estimation takes place after the non-parametric
#    estimation has occurred.
# '''

parametric_estimation <- function(lprobs, dist, kmax, i, j, d, degree, states) {
    ###
    #  This function is used in `fit_dsmm()` when estimation = "parametric".
    #
    #  `lprobs` is a vector of probabilities corresponding to the sojourn times
    #       of the sequence.
    #  `dist` is either NA or one of c('unif', 'pois', 'geom', 'nbinom',
    #      'dweibull')
    #  `kmax` is the maximum sojourn time of the sequence.
    #  `i`, `j`, `d` are values corresponding to the current u, v and sojourn
    #  time distribution (f_(i/d)).
    #  `degree` is the polynomial degree of the drift.
    ###
    if (dist == 'unif') {
        # par1 <- 44
        # lprobs <- sapply(1:k, function(x) if (x <= par1) 1/par1 else 0)
        # freq <- round(lprobs, 6) * 10**6
        # Original method.
        # theta <- tail(which(lprobs != 0), 1)
        # Old method 1.
        # theta <- which(lprobs != 0)
        # theta <- theta[length(theta)]
        # abovezeroprobs <- length(which(lprobs > 0))
        # diffs <- sapply(seq_along(lprobs), function(i) abs(lprobs[i] - lprobs))
        # theta <- median(apply(diffs, c(2),
        #              function(diff) length(which(diff <= 1/abovezeroprobs))))
        # New method 2 - similar consecutive values.
        # ratio <- c(lprobs[-1], NA) / lprobs
        # similar_consec_ratios <- which(ratio <= 2 & ratio >= 0.5)
        # theta <- max(sapply(split(similar_consec_ratios,
        #                       cumsum(c(1, diff(similar_consec_ratios) != 1))),
        #                 length)) + 1
        # New method 3 - similar values regardless of position.
        pos_prob <- lprobs[which(lprobs > 0)]
        all_ratios <- lapply(pos_prob, function(prob) prob / pos_prob)
        valid_ratios <- lapply(all_ratios, function(ratios)
            which( 0.5 <= ratios & ratios <= 2))
        similar_count <- sapply(valid_ratios, length)
        # First approach - Just get the first maximum value...
        # whichmax <- which.max(similar_count)
        # theta <- round(1/ mean(pos_prob[unlist(valid_ratios[whichmax])]))
        # Second approach - Choose among the maximum values the one that
        # returns a sum closer to 1.
        # maxcount <- max(similar_count)
        # theta <- maxcount
        # which_are_max <- which(similar_count == maxcount)
        # mean_thetas <- round(sapply(
        #     which_are_max,
        #     function(max_i) 1 / mean(pos_prob[valid_ratios[[max_i]]])
        # ))
        # diffs <- mean_thetas - maxcount
        # theta <- mean_thetas[which.min(diffs)]
        # Third approach - simply return the number of similar values.
        maxcount <- max(similar_count)
        theta <- maxcount
        return(c(theta, NA))
    } else if (dist == 'pois') {
        # kmax <- 10
        # par1 <- 4 # vector of nonnegative means
        # lprobs <- dpois(0:(kmax-1), lambda = par1)
        # freq <- round(lprobs, 6) * 10**6
        # thetaold <- sum(0:(kmax - 1) * freq) / sum(freq)
        theta <- sum(0:(kmax - 1) * lprobs)
        # theta; thetaold
        return(c(theta, NA))
    } else if (dist == 'geom') {
        # kmax <- 10
        # par1 <- 0.35234  # prob. of success for each trial. (bigger than 0)
        # lprobs <- dgeom(0:(kmax - 1), prob = par1)
        # freq <- round(lprobs, 6) * 10**6
        # thetaold <- sum(freq) / sum(1:kmax * freq)
        theta <- 1 / sum(1:kmax * lprobs) # We cannot divide by 0, thus 1:kmax.
        return(c(theta, NA))
    } else if (dist == 'nbinom') {
        # kmax <- 100
        # par1 <- 10 # number of successful trials - or dispersion parameter.
        # par2 <- 0.7  # probability of success in each trial > 0
        # lprobs <- dnbinom(0:(kmax - 1), size = par1, prob = par2)
        zerokmax <- 0:(kmax - 1)
        expectation <- sum(zerokmax * lprobs)
        expectationx2 <- sum((zerokmax**2) * lprobs)
        expectationsquared <- expectation**2
        variance <- expectationx2 - expectationsquared
        if (expectation >= variance) {
            stop("The negative binomial distribution is not appropriate ",
                 "for modeling the conditional sojourn time distribution ",
                 "associated to the current state u = ", states[i],
                 ", the next state v = ", states[j], " and the distribution ",
                 names_i_d(degree, 'f')[d],
                 ", because variance = ", variance,
                 " >= ", expectation, " = expectation.")
        }
        phat <- expectation / variance
        alphahat <- expectationsquared / (variance - expectation)
        # phat; alphahat
        return(c(alphahat, phat))
        # old
        # req <- round(lprobs, 20) * 10**20
        # xbar <- sum(0:(kmax - 1) * freq) / sum(freq)
        # s2 <- (1 / (sum(freq) - 1)) * sum(freq * (0:(kmax - 1) - xbar)** 2)
        # if (xbar >= s2) {
            # stop("The negative binomial distribution is not appropriate ",
            #     "for modeling the conditional sojourn time distribution ",
            #     "associated to the current state u = ", states[i], ", the next ",
            #     "state v = ", states[j], " and the distribution ",
            #     names_i_d(degree, 'f')[d], ", because variance = ", s2,
            #     " < expectation = ", xbar, ".")
        # }
        # alphahat <- xbar**2 / (s2 - xbar)
        # phat <- xbar / s2
        # theta0 <- c(alphahat, phat)
        # # Constraints about the values of the parameters.
        # # alpha, p > 0.
        # u0 <- diag(x = 1, nrow = 2)
        # c0 <- c(0, 0)
        # # p < 1
        # u1 <- rbind(c(0,-1))
        # c1 <- -1
        # loglik_nbinom <- function(par) {
        #     mask <- which(freq != 0)
        #     kmask <- (0:(kmax - 1))[mask]
        #     fk <- rep.int(x = 0, times = kmax)
        #     fk[mask] <- dnbinom(x = kmask, size = par[1], prob = par[2],
        #                         log = TRUE)
        #     ans <- (-1)*(sum(freq * fk))
        #     ans
        # }
        # mle <- constrOptim(
        #     theta = theta0,
        #     f = loglik_nbinom,
        #     ui = rbind(u0, u1),
        #     ci = c(c0, c1),
        #     method = "Nelder-Mead"
        # )
        # theta <- mle$par
        # return(theta)
    } else if (dist == 'dweibull') {
        if (length(unique(lprobs)) == 1) {
            stop("The Discrete Weibull is not appropriate for modeling",
                 " the conditional sojourn time distribution",
                 " describing the previous state u = ", states[i],
                 ", the next state v = ", states[j],
                 " for the sojourn time distribution ",
                 names_i_d(degree, 'f')[d], ", since the estimation of",
                 " the second parameter beta is not possible.\n",
                 "This happens because we only have 1 non-negative value",
                 " and the estimation of beta requires at least 2 non-negative",
                 " values, which makes the estimation of beta impossible in",
                 " this case")
        }
        qhat <- 1 - lprobs[1]
        cumsumf <- 1 - cumsum(lprobs)
        beta_i <- sapply(2:kmax, function(i)
            log(log(cumsumf[i], base = qhat), base = i))
        betahat <- mean(beta_i[is.finite(beta_i)])
        if (is.nan(betahat)) {
            stop("The Discrete Weibull is not appropriate for modeling",
                 " the conditional sojourn time distribution",
                 " describing the previous state u = ", states[i],
                 ", the next state v = ", states[j],
                 " for the sojourn time distribution ",
                 names_i_d(degree, 'f')[d], ",  the estimation of",
                 " the second parameter beta is not possible.\n ",
                 "This behaviour is perhaps accounted to the fact that",
                 " beta > 1, which is generally hard to estimate.")
        }
        return(c(qhat, betahat))
        # freq <- round(lprobs, 10) * 10**10
        # answer <- tryCatch(
        #     expr = {
        #         # We use suppress warnings because the only warning
        #         # is regarding the times = freq[x] argument (?)
        #         suppressWarnings(
        #             DiscreteWeibull::estdweibull(
        #                 x = unlist(sapply(1:kmax, function(x) rep(x, freq[x]))),
        #                 method = 'ML', zero = FALSE
        #             )
        #         )
        #     }, message = function(m) {
        #         paste0("For the conditional sojourn time distribution",
        #                "describing the previous state u = ", states[i],
        #                ", the next state v = ", states[j],
        #                " and for the sojourn time distribution ",
        #                names_i_d(degree, 'f')[d],
        #                ", function `DiscreteWeibull::estdweibull()`",
        #                "gave the following message:\n\n",
        #                m)
        #     }
        # )
        # if (is.character(answer)) {
        #     error <- paste0(strsplit(answer, ':', fixed = TRUE)[[1]][-4],
        #                     collapse = ":")
        #     stop(error)
        # } else {
        #     theta0 <- answer
        # }
        #
        # # Constraints about the values of the parameters;
        # # q, beta > 0
        # u0 <- diag(x = 1, nrow = 2)
        # c0 <- c(0, 0)
        # # q < 1
        # u1 <- rbind(c(-1,0))
        # c1 <- -1
        # loglik_dweibull <- function(par) {
        #     mask <- which(freq != 0)
        #     kmask <- (1:kmax)[mask]
        #     fk <- rep.int(x = 0, times = kmax)
        #     fk[mask] <- log(
        #         DiscreteWeibull::ddweibull(x = kmask, q = par[1],
        #                                    beta = par[2], zero = FALSE)
        #     )
        #     ans <- (-1)*(sum(freq * fk))
        #     ans
        # }
        # mle <- constrOptim(
        #     theta = theta0,
        #     f = loglik_dweibull,
        #     ui = rbind(u0, u1),
        #     ci = c(c0, c1),
        #     method = "Nelder-Mead"
        # )
        # theta <- mle$par
    }
}

