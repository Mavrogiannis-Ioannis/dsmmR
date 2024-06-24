---
title: NEWS
editor_options: 
  markdown: 
    wrap: 72
---

# dsmmR 1.0.4

## DESCRIPTION

- `fit_dsmm` now has a new attribute, `multi_estimation`, that enables the 
  estimation of a drifting semi-Markov model using multiple sequences. There 
  are two possible options: `avg_model` and `count_sum`.
     - `avg_model` averages the `q_i` received from multiple sequences 
     - `count_sum` adds the counts of the states for each sequence (of equal size)
       and then computes the `q_i`. 
       
- Further changes were made for clarity, to the naming conventions of
  `simulate.dsmm()`. Specifically, the previous variable `seq_length` 
  is now named `max_seq_length`.

# dsmmR 1.0.3

## DESCRIPTION

- Updated the `Depends` section. Now we impose the requirement for R >= 3.5.0,
  in order to make proper use of the `isTRUE()` and `isFALSE()` functions in 
  the `is_logical()` function defined in `utils.R`. These functions will remain
  for their clarity.
  

## Documentation

- Updated the `fit_dsmm()` and `simulate.dsmm()` functions to properly explain the 
  difference between the given/simulated sequence of states and the embedded Markov
  chain. Also, the function `base::rle()` is mentioned for clarity.


# dsmmR 1.0.2

## Documentation

- Updated the `README` and `DESCRIPTION` files with an acknowledgement section.


# dsmmR 1.0.1

## Minor Improvements

-   Added a `NEWS.md` file to track changes to the package.

-   Now the `fit_dsmm()` function has a default value for the `states` attribute,
    being the sorted unique values of the `sequence` character vector attribute.
    

## Bug fixes

-   Fixed a case where `simulate.dsmm()` sometimes did not function as expected
    when `nsim = 1`.
    -   Now it is possible to specify `nsim = 0`, so that the simulated
        sequence will only include the initial state and its
        corresponding sojourn time, e.g. "a", "a", "a".
        
        By giving `nsim = 1` , a single simulation will be made from the
        drifting semi-Markov kernel, returning for example "a", "a",
        "a", "c".

## Documentation

-   Updated the documentation for `simulate.dsmm()`, with accordance to
    the changes made.

-   Updated the `README` file.

    -   Added high-level documentation of the package.
    
    -   Added installation instructions with access to
        the development version of the package through github.

-   Updated the documentation for `dsmmR-package`.

    -   Added a "Community Guidelines" section, so that users can report
        errors or mistakes and contribute directly to the software
        through the newly-established open-source github page at
        <https://github.com/Mavrogiannis-Ioannis/dsmmR>.

    -   Added a "Notes" section, specifying that automated tests are in
        place in order to aid the user with any false input made and,
        furthermore, to ensure that the functions used return the
        expected output.
