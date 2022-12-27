---
title: 'dsmmR: Estimation and Simulation of Drifting Semi-Markov Models'
tags:
	- R
	- semi-Markov Models
	- drifting Markov Models
	- Inhomogeneous Markov chains 
	- DNA sequences
authors:
	- name: Vlad Stefan Barbu
    orcid: 0000-0002-0840-016X
    equal-contrib: true
    affiliation: 1 
	- name: Nicolas Vergne
    equal-contrib: true 
    affiliation: 1
	- name: Ioannis Mavrogiannis
    equal-contrib: true
    affiliation: 1
    affiliations:
	- name: Laboratoire de Mathématiques Raphaël Salem, Université de Rouen Normandie, France
	  index: 1
date: 4 December 2022
bibliography: dsmm_refs.bib
---

# Summary
This is an R package that allows the user to estimate, simulate and define 
different Drifting Semi-Markov model specifications. 
They are used as an alternative against the usual approach with Markov models 
for the Markov models are a usual approach in the modeling of many phenomenons
with a finite state space under discrete-time. However, this imposes the assumption that that the
sequence is homogeneous with respect to time and furthermore that the distribution of the sojourn 
time is Geometric. This is not always true in practice when modeling, for example, DNA sequences.
Drifting Semi-Markov models are build upon Semi-Markov models and Drifting Markov models, allowing
the choice of arbitrary distributions for the sojourn times of the states and at the same time
describe the non-homogeneity through a smooth, known shape that is gradually evolving. In this work,
this gradual evolution is modeled under a polynomial function, which allows the kernel of the Drifting
Semi-Markov models to vary for every instance of the sequence. The degree of the polynomial function 
and the 

The estimation of the Drifting Semi-Markov kernel is non-parametric and three different models can 
be chosen, with the freedom to increase the degree of the polynomial function.

For the parametric Drifting Semi-Markov models specification several discrete sojourn 
time distributions are considered for the sojourn times: Uniform, Geometric, Poisson, 
Discrete Weibull of type 1 and Negative Binomial. The non-parametric specification
makes no assumptions about the shape of the sojourn time distributions. 

The estimation is non-parametric.

Up to three possible model types are allowed for the Drifting Semi-Markov models 
specification and estimation, that concerns whether we allow the Semi-Markov 
kernels to drift with regards to the transition matrix, the sojourn time 
distributions or with both of them together.

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package `astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
by students and experts alike.

# Acknowledgements 


# References
