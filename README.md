# SA-SALSA extension for OpenFOAM
Implementation of the SALSA extension for the Spalart-Allmaras turbulence model in OpenFOAM.

## TODO

- validate implementation
- check why OF throws *Duplicate entry* errors on start-up, but runs normally afterward
- documentation

## References
- template / directory structure taken from https://github.com/TUFRG/SAH-RANS-OF
- T. Rung, U. Bunge, M. Schatz, and F. Thiele, *Restatement of the Spalart–Allmaras Eddy-Viscosity Model 
in Strain-Adaptive Formulation*, AIAA Journal, Vol. 41, no. 7, May 2012, [https://doi.org/10.2514/2.2089](https://arc.aiaa.org/doi/10.2514/2.2089)
- D.-M. Zimmermann, R. Mayer, T. Lutz, and E. Krämer, *Impact of model parameters of SALSA turbulence model
on transonic buffet prediction,* AIAA Journal, Vol. 56, no. 2, pp. 874–877, December 2017,
[https://arc.aiaa.org/doi/10.2514/1.J056193](https://arc.aiaa.org/doi/10.2514/1.J056193)