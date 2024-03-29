% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic_functions.R
\name{simulate.dsmm}
\alias{simulate.dsmm}
\alias{dsmm_simulate}
\title{Simulate a sequence given a drifting semi-Markov kernel.}
\usage{
\method{simulate}{dsmm}(
  object,
  nsim = NULL,
  seed = NULL,
  max_seq_length = NULL,
  klim = 100,
  ...
)
}
\arguments{
\item{object}{An object of S3 class \code{dsmm},
\code{dsmm_fit_nonparametric}, \code{dsmm_nonparametric},
\code{dsmm_fit_parametric} or \code{dsmm_parametric}.}

\item{nsim}{Optional. An integer specifying the number of simulations to be made
from the drifting semi-Markov kernel. The maximum value of \code{nsim} is the
model size which is specified in \code{obj}. This is also the default value.
We define a special case for \code{nsim = 0}, where only the initial distribution
is considered and only the simulation of its sojourn time will be made, without
the next state.}

\item{seed}{Optional. An integer specifying the initialization of the random
number generator.}

\item{max_seq_length}{Optional. A positive integer that will ensure the simulated
sequence will not have a \emph{maximum total length} greater than
\code{max_seq_length} (however, it is possible for the total length to be
\emph{less} than \code{max_seq_length}).}

\item{klim}{Optional. Positive integer. Passed down to \code{get_kernel}
for the parametric object, with class \code{dsmm_parametric}.
Default value is \eqn{100}.}

\item{...}{Optional. Attributes passed down from the \code{simulate} method.}
}
\value{
Returns the simulated sequence for the given drifting
 semi-Markov model. It is a character vector based on \code{nsim} simulations,
 with a maximum length of \code{max_seq_length}.

 This sequence is not to be confused with the embedded Markov chain. The user
 can apply the \code{base::rle()} function on this simulated sequence, if he wishes
 to obtain the corresponding embedded Markov chain and the sojourn times.
}
\description{
Generic function that simulates a sequence under the rule of a
drifting semi-Markov kernel. The number of simulated states is \code{nsim},
while the kernel is retrieved from the object \code{obj} via inheritance
from the S3 class \code{dsmm}.
}
\examples{
# Setup.
sequence <- create_sequence("DNA", len = 1000)
states <- sort(unique(sequence))
d <- 1
obj_model_3 <- fit_dsmm(sequence = sequence,
                        states = states,
                        degree = d,
                        f_is_drifting = TRUE,
                        p_is_drifting = FALSE)

# Using the method `simulate.dsmm()`.
simulated_seq <- simulate(obj_model_3, seed = 1)
short_sim <- simulate(obj = obj_model_3, nsim = 10, seed = 1)
cut_sim <- simulate(obj = obj_model_3, max_seq_length = 10, seed = 1)
str(simulated_seq)
str(short_sim)
str(cut_sim)

# To obtain the embedded Markov chain (EMC) and the corresponding sojourn times
# of any simulated sequence, we can simply use the `base::rle()` function.

sim_seq_emc <- base::rle(cut_sim)$values # embedded Markov chain
sim_seq_sojourn_times <- base::rle(cut_sim)$lengths # sojourn times
cat("Start of the simulated sequence: ", head(cut_sim),
    "...\nThe embedded Markov chain:       ", head(sim_seq_emc),
    "...\nThe sojourn times:               ", head(sim_seq_sojourn_times), "...")

}
\seealso{
About random number generation in R: \code{\link[base:RNG]{RNG}}.

Fitting a model through a sequence from this function: \link{fit_dsmm}.

For the theoretical background of drifting semi-Markov models: \link{dsmmR}.

For obtaining the lengths and values of equals values in a vector:
\code{\link[base:rle]{rle}}.
}
