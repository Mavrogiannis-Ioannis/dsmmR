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
    affiliation: 1
    corresponding: true
  - name: Nicolas Vergne
    affiliation: 1
affiliations:
  - name: Univ Rouen Normandie, CNRS, LMRS UMR 6085, F-76000 Rouen, France
    index: 1
date: 28 December 2022
bibliography: paper.bib
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

Markov models are a common approach in the researcher's toolbox in order to deal with the modeling of many real-life phenomenons presented as a sequence under discrete time. However, this imposes the assumption that the sequence under inspection is homogeneous with respect to time and furthermore that the sojourn times follow the Geometric distribution. This is not always true in practice when modeling, for example, DNA sequences. 

Drifting semi-Markov models (DSMM) formulate the combination of semi-Markov models, which allow arbitrary distributions for the sojourn times, together with drifting Markov models, which describe the non-homogeneity of a sequence through a smooth, known shape that is gradually evolving, expressed through a linear (polynomial) function. As a result, DSMM are best suited to capture the non-homogeneity that occurs from the intrinsic evolution of the system or from the interactions between the system and the environment. For a detailed introduction to semi-Markov models see [@barbu_limnios]. Drifting Markov models were first introduced in [@drift_polynomial].

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


# Acknowledgements

We acknowledge the DATALAB Project (financed by the European Union with the European Regional Development fund (ERDF) and by the Normandy Region) and the HSMM-INCA Project (financed by the French Agence Nationale de la Recherche (ANR) under grant ANR-21-CE40-0005).


# References
