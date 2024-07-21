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
  - name: Laboratoire de Mécanique, Modélisation & Procédés Propres (M2P2), UMR 7340, Aix Marseille Université, CNRS, Centrale Méditerranée, France 
    index: 2
date: 21 July 2024
bibliography: paper.bib
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

Markov models are a common approach in the researcher's toolbox for aiding with the modeling of many real-life phenomena, often presented as a sequence of states in discrete time, known as a Markov chain. However, this assumes that the sequence under inspection is time-homogeneous and that the sojourn times always follow a geometric distribution. Notice that this is not always true in practice, for example when modeling DNA sequences.

Drifting semi-Markov models (DSMM) aim to be a discrete-time pairing of semi-Markov models with drifting Markov models. Semi-Markov models allow an arbitrary choice for the distribution of the sojourn times. The drifting Markov models describe the non-homogeneity of a sequence through a smooth, known shape that is gradually evolving, expressed through a polynomial function. As a result, DSMM are best suited to capture non-homogeneities which occurs from the intrinsic evolution of the system or from the interactions between the system and the environment. For a detailed introduction to semi-Markov models see @barbu_limnios. Drifting Markov models were first introduced in @drift_polynomial.

`dsmmR` is an R package which allows the user to perform parametric and non-parametric estimation and simulation of drifting semi-Markov processes. The user can also define their own parametric and non-parametric DSMM specifications, allowing for a necessary degree of freedom when dealing with a research question. Furthermore, three different types of DSMM are considered. These three models differ in the way they characterize the drifting semi-Markov kernel. Specifically, the first model allows both the transition matrix and the sojourn time distribution to vary (i.e. to "drift"), the second model allows only the transition matrix to drift, while the third model allows only the sojourn time distribution to drift.

The main functions of `dsmmR` are the following:

-   `fit_dsmm()`          : estimates a DSMM given a sequence.
-   `parametric_dsmm()`   : defines a parametric DSMM.
-   `nonparametric_dsmm()`: defines a non-parametric DSMM.
-   `simulate()`          : simulates a sequence given a `dsmm` object.
-   `get_kernel()`        : computes the drifting semi-Markov kernel given a `dsmm` object.

The estimation of the DSMM is parametric or non-parametric and can be defined through the function `fit_dsmm()`, which returns an object of the S3 class `(fit_dsmm_parametric`, `dsmm)` or `(fit_dsmm_nonparametric`, `dsmm)`. In the parametric estimation case, several sojourn time distributions are considered for the sojourn times: Uniform, Geometric, Poisson, Discrete Weibull (of type 1) and Negative Binomial. The parametric DSMM specification can be defined through the function `parametric_dsmm()`, which returns an object of the S3 class `(dsmm_parametric`, `dsmm)`. The non-parametric DSMM specification can be defined through the function `nonparametric_dsmm()`, which returns an object of the S3 class `(dsmm_nonparametric`, `dsmm)`. It allows for the sojourn time distributions to be of an arbitrary shape.

The `dsmm` class acts like a wrapper class, enabling the handling of all above objects, which encompass the three models and parametric or non-parametric cases for any degree for:

-   Simulating a sequence of states under a drifting semi-Markov kernel through the S3 method `simulate.dsmm()`.

-   Calculating the drifting semi-Markov kernel through the generic function `get_kernel()`.

# Statement of need

The present R package consists of a novel approach to model the gradually evolving non-homogeneity of a discrete-time Markov chain through a polynomial function, while it simultaneously allows for the flexibility of choosing arbitrary sojourn time distributions. Drifting semi-Markov processes represent a versatile modeling tool, which can be readily applied in many fields of science like bio-informatics, engineering, finance and also potentially in reliability and survival analysis settings. In the literature, existing R packages deal with semi-Markov models, hidden semi-Markov models and drifting Markov models.

Such R packages specializing in semi-Markov models include the `SemiMarkov` [@SemiMarkov_package] package, which performs maximum likelihood estimation for parametric continuous-time semi-Markov processes, while they allow a choice of distribution for the sojourn times between the Exponential, the Weibull or the exponentiated Weibull. Also associated hazard rates are computed, and covariates can also be taken into account through the Cox proportional hazard model. The R package `smmR` [@smmR_package] deals with discrete-time multi-state semi-Markov models, performing parametric and non-parametric estimation and simulation, with k-th order Markov chains also considered. For the parametric estimation the following sojourn times were considered: Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial. One or more sample paths were considered, with or without censoring at the beginning or the end of the sample paths. There are also a few R packages also dedicated to hidden semi-Markov models, where methods for estimation and prediction are implemented, such as the `HMM` R package [@HMM_package], the `HiddenMarkov` R package [@HiddenMarkov_package] and the `mhsmm` R package [@mhsmm_package]. 

Furthermore, for drifting Markov models the R package `drimmR` [@drimmR_package] has been developed, performing non-parametric estimation and simulation, also allowing for the exact computation of the associated reliability measures such as reliability, availability, maintainability and failure rates. Also several statistical frameworks were considered, allowing for estimation from multiple sequences, complete or incomplete, with the same or unequal model size. In this package, a polynomial function was also used to capture the non-homogeneity of the sequences considered.


# Acknowledgements

We acknowledge the DATALAB Project (financed by the European Union with the European Regional Development fund (ERDF) and by the Normandy Region) and the HSMM-INCA Project (financed by the French Agence Nationale de la Recherche (ANR) under grant ANR-21-CE40-0005).


# References
