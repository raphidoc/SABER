# S.A.B.E.R
## Semi Analytical Bayesian Error Retrieval

A radiative transfer model to retrieve inherent optical proprieties, water depth and bottom reflectance from remote sensing reflectance.
Implement Markov Chain Monte Carlo

Creator and original developer: Soham Mukherjee
Packaging: Raphael Mabit

## Questions

* Where does the pure water IOPs comes from ? a = Pope 97, b = Morel 74
* In pure water IOP Why different pure water bb values for case one or two waters ? pure water IOP does not depends on water type.
* Where does the parametric formula to retrieve spectral slope of CDOM + NAP comes from, QAA ? QAA v5 or 6
* Use `saber_forward_parametric_conc_wise.R` as authoritative final function.
* The inversion of fraction rb_class is not constrained to unity ? Seems like a difficult problem ...
* What about the geometry in inversion ? With whihch value do we parametrize forward models ?
* Is it really rrs_bottom [sr^-1], is it not just Rb [dimensionless] ? Check the SVC method !
* What's the logic between negative and positive log-ll? minimizing the error vs maximizing the likelihood?
