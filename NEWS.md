---
editor_options: 
  markdown: 
    wrap: 72
---

# dsmmR 1.0.1

## Minor Improvements

-   Added a `NEWS.md` file to track changes to the package.

## Bug fixes

-   Fixed a case where `simulate.dsmm()` did not function as expected
    when `nsim = 1`.
    -   Now it is possible to specify `nsim = 0`, so that the simulated
        sequence will only include the initial state and its
        corresponding sojourn time, e.g. "a", "a", "a".\
        By giving `nsim = 1` , a single simulation will be made from the
        drifting semi-Markov kernel, returning for example "a", "a",
        "a", "c".

## Documentation

-   Updated the documentation for `simulate.dsmm()`, with accordance to
    the changes made.

-   Updated the `README` in the following manner:

    -   Added high-level documentation of the package alongside with
        proper installation instructions, as well as giving access to
        the development version of the package through github.

    -   Added a "Community Guidelines" section, so that users know where
        to report errors or mistakes and allowing them to contribute
        directly to the software through the newly-established open
        source github page.

    -   Added a "Notes" section, specifying that automated tests are in
        place in order to aid the user with any false input made and,
        furthermore, to ensure that the functions used return the
        expected output.
