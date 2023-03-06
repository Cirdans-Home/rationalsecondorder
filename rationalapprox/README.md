# Rational Approximation

This folder contains the codes related to the different rational 
approximations of the sinc function and the generation of different 
types of poles.

- `genlagpol.m` Laguerre based poles,
- `genlagpolsym.m` Doubling of the Laguerre based poles,
- `genlagexppol.m` Laguerre based poles from the exponential expansion,
- `genpadeexppol.m` Padé poles for the exp(-z) function,
- `cf.m` Poles and residual of the Carathéodory-Fejér rational approximation to the exponential,
- `sincfourier.m` Computation of `sinc(A)v` via inversion of the Fourier transform (exponential sums),
- `sincm.m` Rational Krylov method for the `sinc(A)v` computation,
- `testingwithtimings.m`, `testofpoles.m`, `plotresults.m` execute the tests in the paper.
