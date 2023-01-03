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
  - name: Laboratoire de Mathématiques Raphaël Salem, Université de Rouen Normandie, France
    index: 1
date: 28 December 2022
bibliography: paper.bib
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

## General Introduction for a non-specialist audience.

Markov models are a common approach in the researcher's toolbox in order to deal with the modeling of many real-life phenomenons presented as a sequence under discrete time. However, this imposes the assumption that the sequence under inspection is homogeneous with respect to time and furthermore that the sojourn times follow the geometric distribution. This is not always true in practice when modeling, for example, DNA sequences. 

### Purpose

Drifting semi-Markov models (DSMM) formulate the combination of semi-Markov models (see @barbu_limnios for a detailed introduction), which allow arbitrary distributions for the sojourn times, together with the drifting Markov models (introduced first in @drift_polynomial), which describe the non-homogeneities of a sequence through a smooth, known shape that is gradually evolving, expressed through a linear (polynomial) function.

As a result, DSMM are best suited to capture non-homogeneities that occur from the intrinsic evolution of the system or from the interactions between the system and the environment.

`dsmmR` is an R package which allows the user to perform parametric and non-parametric estimation and simulation of drifting semi-Markov processes. The user can also define their own parametric and non-parametric DSMM specifications, allowing for a necessary degree of freedom when dealing with a research question.


These models differ in the number of transition matrices and sojourn time distributions used for the computation of a number of semi-Markov kernels, which in turn characterize the Drifting semi-Markov kernel. 

Furthermore, three different types of DSMM are considered. These three models differ in the way they characterize the drifting semi-Markov kernel. Specifically, the first model allows both the transition matrix and the sojourn time distribution to vary (i.e. to "drift"), the second model allows only the transition matrix to drift while the third model allows only the sojourn time distribution to drift.

### High-level documentation

The main functions of `dsmmR` are the following:

-   `fit_dsmm()`: estimates a DSMM given a sequence.
-   `parametric_dsmm()`: defines a parametric DSMM.
-   `nonparametric_dsmm()`: defines a non-parametric DSMM.
-   `simulate()`: simulates a sequence.
-   `get_kernel()`: returns the drifting semi-Markov kernel.

The estimation of the DSM kernel is parametric or non-parametric and can be defined through the function `fit_dsmm()`, which returns an object of the S3 class (`fit_dsmm_parametric`, `dsmm`) or (`fit_dsmm_nonparametric`, `dsmm`). In the parametric estimation case, several discrete sojourn time distributions are considered for the sojourn times: uniform, geometric, poisson, discrete weibull (type 1) and negative binomial.



The parametric DSMM specification can be defined through the function `parametric_dsmm()`, which returns an object of the S3 class (`dsmm_parametric`, `dsmm`). In the parametric definition case the same discrete sojourn time distributions. 

The non-parametric DSMM specification can be defined through the function `nonparametric_dsmm()`, which returns an object of the S3 class (`dsmm_nonparametric`, `dsmm`). It allows for the sojourn time distributions to be of an arbitrary shape.


## Perhaps not necessary; we speak of the 3 models above.
For all three of the above functions, the user can define if they want the transition matrix and the sojourn times to be drifting (model 1), if they want only the transition matrix to drift (model 2) or if they want only the sojourn time distribution to drift (model 3), resulting in three possible modeling options. Also, a choice can be made with regards to the degree of the drift. With the DSMM definitions, one can also define the model size (length of the EMC).

Based on any of these classes, the following methods (generic functions)
are available:

-   Simulate a sequence of states under a drifting semi-Markov kernel
    through the function `simulate.dsmm()`.

-   Obtain the drifting semi-Markov kernel through the generic function
    `get_kernel()`.

Thus, the class `dsmm` acts like a wrapper class for drifting
semi-Markov model specifications.


# Statement of need

The semi-Markov processes represent a versatile tool that is applied in many fields of science like reliability, survival analysis, bioinformatics, engineering, finance and more.

It is necessary to enable science a new path in order to deal with the non-homogeneities of sequences that might be present in a research question.

Few R packages have been developed to handle semi-Markov models or hidden semi-Markov models. For semi-Markov models we have the recent semi-Markov R package (Listwon & Saint-Pierre, 2015) that performs maximum likelihood estimation for parametric continuous-time semi-Markov processes, where the distribution can be chosen between Exponential, Weibull or exponentiated Weibull. That package computes associated hazard rates; covariates can also be taken into account through the Cox proportional hazard model. Few R packages are also dedicated to hidden semi-Markov models, implementing estimation and prediction methods. Among them, we can cite the hsmm R package (Bulla et al., 2010) and the mhsmm R package (O'Connell et al., 2011). The package SMM (Barbu et al., 2018) deals with discrete-time multi-state semi-Markov models but does not.

# Estimation

# Acknowledgements

This was part of the DataLab Normandie project.

We acknowledge the project AStERiCs Apprentissage Statistique à l'Echelle pour la Représen- tation et la Classification non-supervisées (RIN project funded by the Normandy Region), DAISI on Biomedical Data Classification (co-financed by the European Union with the European Regional Development Fund (ERDF) and by the Normandy Region). We also acknowledge Mathilde Sautreuil, Caroline Bérard and Dominique Cellier for the help they provided in creating the first working package SMM (Barbu et al., 2018).

# References
