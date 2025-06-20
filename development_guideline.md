# External dependencies (saber-lib and cmdstan)

The SABER R package is a high level wrapper around saber-lib and cmdstan.
To work properly, it need a local installation of both those software.

saber-lib contains all the bio-optical and radiative transfer models.
cmdstan is used to compile a STAN model for Markov Chain Monte Carlo inversion.

# autoconf

*Check and document the usage of autoconf in this package*

Modify file configure.ac and Makevars.in then run:

```
autoreconf -fi -I m4
./configure
```
