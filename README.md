# S.A.B.E.R
## Semi Analytical Bayesian Error Retrieval

A radiative transfer model to retrieve inherent optical proprieties, water depth and bottom reflectance from remote sensing reflectance.
Implement Markov Chain Monte Carlo

Creator and original developer: Soham Mukherjee
Packaging: Raphael Mabit

## Outline of the code structure

This package follow the recommendations of https://r-pkgs.org/ and the tidyverse style guide https://style.tidyverse.org/.

The code is written with a functional approach. The "business" logic
(the low level functions) are written in file names starting with (`fct_*`),
more generic function are written under (`utils_*`) files.

The central piece of code that stitches together the forward models with the objective function is `objective_factory` high level function.
It allows to easily combine any forward model with any objective functions, even those that a user might add.

For users to add their own forward models and objectives function, we adopted the use of registries.
A registry is an environment in which specific function are stored and can
be retrieved by name in the higher level function arguments.
The three registry currently in use are `.input_preparer_registry`, `.forward_model_registry`, `.objective_function_registry`.

## Questions

* Where does the pure water IOPs comes from ? a = Pope 97, b = Morel 74
* In pure water IOP Why different pure water bb values for case one or two waters ? pure water IOP does not depends on water type.
* Where does the parametric formula to retrieve spectral slope of CDOM + NAP comes from, QAA ? QAA v5 or 6
* Use `saber_forward_parametric_conc_wise.R` as authoritative final function.
* The inversion of fraction rb_class is not constrained to unity ? Seems like a difficult problem ...
* What about the geometry in inversion ? With whihch value do we parametrize forward models ?
* Is it really rrs_bottom [sr^-1], is it not just Rb [dimensionless] ? Check the SVC method !
* What's the logic between negative and positive log-ll? minimizing the error vs maximizing the likelihood?

* We should probably not use extrapolation of spectrum. for example the phyto absorption spectrum in `iop_from_oac`
