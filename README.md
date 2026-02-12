# SA-SALSA extension for OpenFOAM
Implementation of the strain-adaptive linear Spalart-Allmaras (SALSA) turbulence model in OpenFOAM.

The repositories of [TUFRG](https://github.com/TUFRG/SAH-RANS-OF) and 
[mAlletto](https://gitlab.com/mAlletto/openfoamtutorials/-/tree/master/SpalartAllmarasRCsend)
were used as template for the structure of this implementation.


## TODO

- validate implementation -> difference $\tilde{C}_{b1}$ and $C_{b1}$ not mentioned in paper
(definition of $\tilde{C}_{b1}$ missing), also in SA $(\tilde{\nu} + \nu) / \sigma$ but in SALSA $\nu + \tilde{\nu} / \sigma$
- documentation of all equations etc.
- check if we can implement this more efficiently

## References
- template / directory structure taken from [TUFRG](https://github.com/TUFRG/SAH-RANS-OF) and 
[mAlletto](https://gitlab.com/mAlletto/openfoamtutorials/-/tree/master/SpalartAllmarasRCsend)
- T. Rung, U. Bunge, M. Schatz, and F. Thiele, *Restatement of the Spalart–Allmaras Eddy-Viscosity Model 
in Strain-Adaptive Formulation*, AIAA Journal, Vol. 41, no. 7, May 2012, [https://doi.org/10.2514/2.2089](https://arc.aiaa.org/doi/10.2514/2.2089)
- D.-M. Zimmermann, R. Mayer, T. Lutz, and E. Krämer, *Impact of model parameters of SALSA turbulence model
on transonic buffet prediction,* AIAA Journal, Vol. 56, no. 2, pp. 874–877, December 2017,
[https://arc.aiaa.org/doi/10.2514/1.J056193](https://arc.aiaa.org/doi/10.2514/1.J056193)