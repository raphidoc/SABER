# S.A.B.E.R
## Semi Analytical Bayesian Error Retrieval

A radiative transfer model to retrieve inherent optical proprieties, water depth and bottom reflectance from remote sensing reflectance.
Implement Markov Chain Monte Carlo

Creator and original developer: Soham Mukherjee
Packaging: Raphael Mabit

## Questions

* Where does the pure water IOPs comes from ?
* Why different pure water bb values for case one or two waters ? pure water IOP does not depends on water type.
* SSR objective function is define in the parameter of inversion but not implemented in the code ?
* Does the function `Saber_forward_fast` actually implement the am03 SA model ?
* Where does the parametric formula to retrieve spectral slope of CDOM + NAP comes from, QAA ?
* What's the difference between `SABER_forward_fast.R` and `saber_forward_parametric_conc_wise.R` ?
* Is Saber_forward the am03 model (is it ok to change the name) ?
* The inversion of fraction rb_class is not constrained to unity ? Seems like a difficult problem ...
* What about the geometry in inversion ? With wihch value do we parametrize forward models ?
* Is it really rrs_bottom [sr^-1], is it not just Rb [dimensionless] ?
* What is the logic of `runMCMC` post processing ? 3 Chains ? Take the mean ?
