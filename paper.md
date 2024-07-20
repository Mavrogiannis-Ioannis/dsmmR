---
title: 'dsmmR: Estimation and Simulation of Drifting Semi\-Markov Models'
tags:
  - R
  - non\-homogeneous Markov chains 
  - semi\-Markov models
  - drifting Markov models
  - DNA sequences
authors:
  - name: Vlad Stefan Barbu
    orcid: 0000-0002-0840-016X
    affiliation: 1 
  - name: Ioannis Mavrogiannis
    affiliation: 2
    corresponding: true
  - name: Nicolas Vergne
    affiliation: 1
affiliations:
  - name: Laboratory of Mathematics Raphaël Salem (LMRS), UMR CNRS 6085, University of Rouen Normandy, France
    index: 1
  - name: Laboratoire de Mécanique, Modélisation & Procédés Propres (M2P2), M2P2 UMR 7340, Aix Marseille Université, CNRS, Centrale Méditerranée, France 
    index: 2
  - name: Univ Rouen Normandie, CNRS, LMRS UMR 6085, F-76000 Rouen, France
    index: 3
date: 20 July 2024
bibliography: paper.bib
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

Markov models are a common approach in the researcher's toolbox in order to deal with the modeling of many real-life phenomenons presented as a sequence under discrete time. However, this imposes the assumption that the sequence under inspection is homogeneous with respect to time and furthermore that the sojourn times follow the Geometric distribution. This is not always true in practice when modeling, for example, DNA sequences.

Drifting semi-Markov models (DSMM) formulate the combination of semi-Markov models, which allow arbitrary distributions for the sojourn times, together with drifting Markov models, which describe the non-homogeneity of a sequence through a smooth, known shape that is gradually evolving, expressed through a linear (polynomial) function. As a result, DSMM are best suited to capture the non-homogeneity that occurs from the intrinsic evolution of the system or from the interactions between the system and the environment. For a detailed introduction to semi-Markov models see @barbu_limnios. Drifting Markov models were first introduced in @drift_polynomial.

`dsmmR` is an R package which allows the user to perform parametric and non-parametric estimation and simulation of drifting semi-Markov processes. The user can also define their own parametric and non-parametric DSMM specifications, allowing for a necessary degree of freedom when dealing with a research question. Furthermore, three different types of DSMM are considered. These three models differ in the way they characterize the drifting semi-Markov kernel. Specifically, the first model allows both the transition matrix and the sojourn time distribution to vary (i.e. to "drift"), the second model allows only the transition matrix to drift while the third model allows only the sojourn time distribution to drift.

The main functions of `dsmmR` are the following:

-   `fit_dsmm()`: estimates a DSMM given a sequence.
-   `parametric_dsmm()`: defines a parametric DSMM.
-   `nonparametric_dsmm()`: defines a non-parametric DSMM.
-   `simulate()`: simulates a sequence.
-   `get_kernel()`: returns the drifting semi-Markov kernel.

The estimation of the DSMM kernel is parametric or non-parametric and can be defined through the function `fit_dsmm()`, which returns an object of the S3 class `(fit_dsmm_parametric`, `dsmm)` or `(fit_dsmm_nonparametric`, `dsmm)`. In the parametric estimation case, several discrete sojourn time distributions are considered for the sojourn times: Uniform, Geometric, Poisson, Discrete Weibull (of type 1) and Negative Binomial. The parametric DSMM specification can be defined through the function `parametric_dsmm()`, which returns an object of the S3 class `(dsmm_parametric`, `dsmm)`. In the parametric definition case, we have the same discrete sojourn time distributions. The non-parametric DSMM specification can be defined through the function `nonparametric_dsmm()`, which returns an object of the S3 class `(dsmm_nonparametric`, `dsmm)`. It allows for the sojourn time distributions to be of an arbitrary shape.

The `dsmm` class acts like a wrapper class, enabling the handling of all three models and parametric or non-parametric cases for any degree when it comes to the following actions:

-   Simulating a sequence of states under a drifting semi-Markov kernel through the S3 method `simulate.dsmm()`.

-   Calculating the drifting semi-Markov kernel through the generic function `get_kernel()`.

# Statement of need

Drifting semi-Markov processes represent a versatile tool that can be applied in many fields of science like reliability, survival analysis, bioinformatics, engineering, finance and more. The present R package consists of a novel approach to tackle gradually evolving non-homogeneity in a polynomial way and also allows for arbitrary sojourn time distributions. Instead, existing R packages deal with semi-Markov models, hidden semi-Markov models and drifting Markov models.

For semi-Markov models we have the R package `SemiMarkov` [@SemiMarkov_package] that performs maximum likelihood estimation for parametric continuous-time semi-Markov processes, where the distribution can be chosen between Exponential, Weibull or exponentiated Weibull. That package computes associated hazard rates; covariates can also be taken into account through the Cox proportional hazard model. The R package `smmR` [@smmR_package] deals with discrete-time multi-state semi-Markov models, performing parametric and non-parametric estimation and simulation, with k-th order Markov chains also considered. For the parametric estimation the following sojourn times were considered: Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial. One or more sample paths were considered, with or without censoring at the beginning or the end of the sample paths. Few R packages are also dedicated to hidden semi-Markov models, implementing estimation and prediction methods. Among them, we can cite the HMM R package [@HMM_package], the `HiddenMarkov` R package [@HiddenMarkov_package] and the `mhsmm` R package [@mhsmm_package]. 

Furthermore, for drifting Markov models the R package `drimmR` [@drimmR_package] was developed, performing non-parametric estimation and simulation, also allowing for the exact computation of the associated reliability measures and for several frameworks with regards to how many samples were used, considering if they were complete or incomplete and if they were of the same length. In this package, a linear (polynomial) function was used to capture the non-homogeneity. 

# Quickstart

```{r, eval = TRUE, include = FALSE}
library(dsmmR)
states <- c("a", "b", "c")
s <- length(states)
degree <- 1
p_dist_1 <- matrix(c(0,   0.4,  0.6,
                     0.5, 0,    0.5,
                     0.3, 0.7,  0   ), ncol = s, byrow = TRUE)
p_dist_2 <- matrix(c(0,   0.55, 0.45,
                     0.25, 0,   0.75,
                     0.5, 0.5,  0   ), ncol = s, byrow = TRUE)
p_dist <- array(c(p_dist_1, p_dist_2), dim = c(s, s, degree + 1))
f_dist_1 <- matrix(c(NA,   "nbinom",   "unif",
                   "geom",  NA,        "pois",
                   "pois", "dweibull",  NA   ), nrow = s, ncol = s, byrow = TRUE)
f_dist_1_pars_1 <- matrix(c(NA,  4,   3,
                            0.7, NA,  5,
                            3,   0.6, NA), nrow = s, ncol = s, byrow = TRUE)
f_dist_1_pars_2 <- matrix(c(NA,  0.5, NA,
                            NA,  NA,  NA,
                            NA,  0.8, NA), nrow = s, ncol = s, byrow = TRUE)
f_dist_2 <- f_dist_1 
f_dist_2_pars_1 <- matrix(c(NA,  3,   5,
                            0.3, NA,  2,
                            5,   0.3, NA), nrow = s, ncol = s, byrow = TRUE)
f_dist_2_pars_2 <- matrix(c(NA,  0.4, NA,
                            NA,  NA,  NA,
                            NA,  0.5, NA), nrow = s, ncol = s, byrow = TRUE)`.

f_dist <- array(c(f_dist_1, f_dist_2), dim = c(s, s, degree + 1))
f_dist_pars <- array(c(f_dist_1_pars_1, f_dist_1_pars_2,
                       f_dist_2_pars_1, f_dist_2_pars_2), 
                     dim = c(s, s, 2, degree + 1))
dsmm_model <- parametric_dsmm(
    model_size = 10000,
    states = states,
    initial_dist = c(0.6, 0.3, 0.1),
    degree = degree,
    p_dist = p_dist,
    f_dist = f_dist,
    f_dist_pars = f_dist_pars,
    p_is_drifting = TRUE,
    f_is_drifting = TRUE
)
sim_seq <- simulate(dsmm_model, klim = 30, seed = 1)
fitted_model <- fit_dsmm(sequence = sim_seq,
                         states = states,
                         degree = degree,
                         f_is_drifting = TRUE,
                         p_is_drifting = TRUE,
                         estimation = 'parametric',
                         f_dist = f_dist)
```


To install the package, simply use either of the commands in the R console:
```{r, eval = FALSE}
install.packages('dsmmR') # official release in CRAN
```
Otherwise, if you wish to download the development version directly from Github: 
```{r, eval = FALSE}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Mavrogiannis-Ioannis/dsmmR")
```

We can then set up a simple case of defining a parametric `dsmm` object,
simulating from it and then estimating from the simulated sequence, 
by first of all loading the package, 
```{r, eval = FALSE}
library(dsmmR)
```
and then defining the states and a degree equal to 1.
```{r, eval = FALSE}
states <- c("a", "b", "c")
s <- length(states)
degree <- 1
```

Since degree is equal to 2, we can then define the 2 drifting transition matrices:
```{r, eval = FALSE}
p_dist_1 <- matrix(c(0,   0.4,  0.6,
                     0.5, 0,    0.5,
                     0.3, 0.7,  0   ), ncol = s, byrow = TRUE)
p_dist_2 <- matrix(c(0,   0.55, 0.45,
                     0.25, 0,   0.75,
                     0.5, 0.5,  0   ), ncol = s, byrow = TRUE)
p_dist <- array(c(p_dist_1, p_dist_2), dim = c(s, s, degree + 1))
```

Let us also consider the case where only the parameters of the distributions 
modeling the sojourn times are drifting across the sequence. Note that some 
distributions like the Negative Binomial and the Discrete Weibull require two 
parameters, which we define in two matrices.
```{r, eval = FALSE}
f_dist_1 <- matrix(c(NA,   "nbinom",   "unif",
                   "geom",  NA,        "pois",
                   "pois", "dweibull",  NA   ), nrow = s, ncol = s, byrow = TRUE)
f_dist_1_pars_1 <- matrix(c(NA,  4,   3,
                            0.7, NA,  5,
                            3,   0.6, NA), nrow = s, ncol = s, byrow = TRUE)
f_dist_1_pars_2 <- matrix(c(NA,  0.5, NA,
                            NA,  NA,  NA,
                            NA,  0.8, NA), nrow = s, ncol = s, byrow = TRUE)
f_dist_2 <- f_dist_1 
f_dist_2_pars_1 <- matrix(c(NA,  3,   5,
                            0.3, NA,  2,
                            5,   0.3, NA), nrow = s, ncol = s, byrow = TRUE)
f_dist_2_pars_2 <- matrix(c(NA,  0.4, NA,
                            NA,  NA,  NA,
                            NA,  0.5, NA), nrow = s, ncol = s, byrow = TRUE)`.

f_dist <- array(c(f_dist_1, f_dist_2), dim = c(s, s, degree + 1))
f_dist_pars <- array(c(f_dist_1_pars_1, f_dist_1_pars_2,
                       f_dist_2_pars_1, f_dist_2_pars_2), 
                     dim = c(s, s, 2, degree + 1))
```

Then, defining a `dsmm_parametric` object is done simply through the function 
`parametric_dsmm()`:
```{r, eval = FALSE}
dsmm_model <- parametric_dsmm(
    model_size = 10000,
    states = states,
    initial_dist = c(0.6, 0.3, 0.1),
    degree = degree,
    p_dist = p_dist,
    f_dist = f_dist,
    f_dist_pars = f_dist_pars,
    p_is_drifting = TRUE,
    f_is_drifting = TRUE
)
```

We can then simulate a sequence from this parametric object like-so:
```{r, eval = FALSE}
sim_seq <- simulate(dsmm_model, klim = 30, seed = 1)
```

To fit this sequence with a drifting semi-Markov model, one can use:
```{r, eval = FALSE}
fitted_model <- fit_dsmm(sequence = sim_seq,
                         states = states,
                         degree = degree,
                         f_is_drifting = TRUE,
                         p_is_drifting = TRUE,
                         estimation = 'parametric',
                         f_dist = f_dist)
```

Finally, the drifting transition matrix is estimated as:
```{r, eval = TRUE}
print(fitted_model$dist$p_drift, digits = 2)
```
and the parameters for the drifting sojourn time distributions are:
```{r, eval = TRUE}
print(fitted_model$dist$f_drift_parameters, digits = 2)
```

For more details about any of these functions, consider viewing the extended
documentation by adding a `?` in front of the function you want to view, 
for example `?fit_dsmm`.


# Acknowledgements

We acknowledge the DATALAB Project (financed by the European Union with the European Regional Development fund (ERDF) and by the Normandy Region) and the HSMM-INCA Project (financed by the French Agence Nationale de la Recherche (ANR) under grant ANR-21-CE40-0005).


# References
