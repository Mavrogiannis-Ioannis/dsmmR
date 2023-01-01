---
title: 'dsmmR: Estimation and Simulation of Drifting Semi\-Markov Models'

tags:
  - R
  - semi\-Markov Models
  - drifting Markov Models
  - non\-homogeneous Markov chains 
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

---

# Summary

Semi\-Markov models are a usual approach in the modelling of many phenomenons with a finite state space under discrete\-time. However, this imposes the assumption that that the sequence is homogeneous with respect to time and furthermore that the distribution of the sojourn times is Geometric. This is not always true in practice when modelling, for example, DNA sequences.

`dsmmR` is an R package that allows the user to estimate, simulate and define parametric and nonparametric drifting semi\-Markov model specifications. 

Drifting semi\-Markov models formulate the combination of semi\-Markov models (see [@barbu_limnios] for an introduction in discrete\-time) together with drifting Markov models (first introduced in [@drift_polynomial]). This allows for the choice of arbitrary distributions for the sojourn times and at the same time it enables the description of the non\-homogeneities present in the sequence through a smooth, gradually evolving shape. 

In this work, this shape is defined through a polynomial function $A_i(t)$ with degree $d$. 
This allows the drifting semi\-Markov kernel kernel $q_{\frac{t}{n}}$ to vary (i.e. to "drift")for every instance $t$ of the embedded Markov chain $\large( J_t \large)_{t \in \{1, \dots, n\}}$ of the sequence. This is achieved through a number of semi-Markov kernels $q_{\frac{i}{d}}$ that are fixed alongside the sequence:

$$q_{\frac{t}{n}} = \sum_t^{+\infty}q_{\frac{i}{d}}A_i(t)$$

### extra info:
The way to conduct a non\-homogeneous fit of the sequence is achieved through the embedded Markov chain (EMC). For a degree equal to $d$, we fix $d + 1$ semi-Markov kernels alongside the sequence, in a uniform way. For example, if $d = 1$, we will have two semi-Markov kernels where we fix one at the beginning of the EMC and one at the end of the EMC. Then, for every instance between these points, a linear functionThe degree also If the degree is equal to 1, the polynomial function will be linear.

it fixes evenly a number of points alongside the sequence, also defines the number of the sojourn time distributions and the transition matrices used, allowing them to “drift” between the same number of points on the sequence.

Two more cases are considered, when only one of the sojourn time distributions or the transition matrices are drifting and the other remains constant.

# Estimation
The estimation of the drifting semi\-Markov kernel is non-parametric and can be defined through the function `fit_dsmm()`, which returns an object of the S3 class (`fit_dsmm`, `dsmm`).

The parametric drifting semi\-Markov model specification can be defined from the user through the function `parametric_dsmm()`, which returns an object of the S3 class (`dsmm_parametric`, `dsmm`). Several discrete sojourn time distributions are considered for the sojourn times: Uniform, Geometric, Poisson, Discrete Weibull of type 1 and Negative Binomial. 

The non-parametric drifting semi\-Markov model specification can be defined from the user through the function `nonparametric_dsmm()`, which returns an object of the S3 class (`dsmm_nonparametric`, `dsmm`). It allows for the sojourn time distributions to be of an arbitrary shape.

For all three of these functions, three possible options for the drift are available, as well as the choice of the degree of the polynomial function of the model.

Based on any of these classes, the following methods are available:

* Simulate a sequence of states under a drifting semi\-Markov kernel through the function `simulate.dsmm()`;

* Obtain the drifting semi\-Markov kernel through the generic function `get_kernel()`.

Thus, the class `dsmm` acts like a wrapper class for drifting semi\-Markov model specifications.


# Summary

The estimation of the drifting semi\-Markov kernel is non-parametric and three different models can be chosen, with the freedom to increase the degree of the polynomial function.

For the parametric drifting semi\-Markov models specification several discrete sojourn time distributions are considered for the sojourn times: Uniform, Geometric, Poisson, Discrete Weibull of type 1 and Negative Binomial. The non-parametric specification makes no assumptions about the shape of the sojourn time distributions.

The estimation is non-parametric.

Up to three possible model types are allowed for the drifting semi\-Markov models specification and estimation, that concerns whether we allow the semi\-Markov kernels to drift with regards to the transition matrix, the sojourn time distributions or with both of them together.

About semi\-Markov models: [@barbu_limnios]
About drifting Markov models: [@drift_polynomial]

The forces on stars, galaxies, and dark matter under external gravitational fields lead to the dynamical evolution of structures in the universe. The orbits of these bodies are therefore key to understanding the formation, history, and future state of galaxies. The field of "galactic dynamics," which aims to model the gravitating components of galaxies to study their structure and evolution, is now well-established, commonly taught, and frequently used in astronomy. Aside from toy problems and demonstrations, the majority of problems require efficient numerical tools, many of which require the same base code (e.g., for performing numerical orbit integration).

# Statement of need

The semi\-Markov processes represent a versatile tool that is applied in many fields of science
like reliability, survival analysis, bioinformatics, engineering, finance, etc. 

It is necessary to enable science a new path in order to deal with the homogeneities
of the sequences that might be present in a research question. 

Few R packages
have been developed to handle semi\-Markov models or hidden semi\-Markov models. For
semi\-Markov models we have the recent semiMarkov R package (Listwon & Saint-Pierre, 2015)
that performs maximum likelihood estimation for parametric continuous-time semi\-Markov
processes, where the distribution can be chosen between Exponential, Weibull or exponentiated
Weibull. That package computes associated hazard rates; covariates can also be taken into
account through the Cox proportional hazard model. Few R packages are also dedicated to
hidden semi\-Markov models, implementing estimation and prediction methods. Among them,
we can cite the hsmm R package (Bulla et al., 2010) and the mhsmm R package (O’Connell et al.,
2011). The package SMM (Barbu et al., 2018) deals with discrete-time multi-state semi\-Markov
models but does not


# Acknowledgements

This was part of the DataLab Normandie project.

We acknowledge the project AStERiCs Apprentissage Statistique à l’Echelle pour la Représen-
tation et la Classification non-supervisées (RIN project funded by the Normandy Region), DAISI
on Biomedical Data Classification (co-financed by the European Union with the European
Regional Development Fund (ERDF) and by the Normandy Region).
We also acknowledge Mathilde Sautreuil, Caroline Bérard and Dominique Cellier for the help
they provided in creating the first working package SMM (Barbu et al., 2018


# References
