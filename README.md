# SA-SALSA extension for OpenFOAM
Implementation of the SALSA extension for the Spalart-Allmaras turbulence model in OpenFOAM.

## TODO

- check implementation of $C_{b1}$ requires adjustment of $C_{w1}$
- validate implementation
- check why OF throws *Duplicate entry* errors on start-up, but runs normally afterward
- documentation

## References
- template / directory structure taken from https://github.com/TUFRG/SAH-RANS-OF
- T. Rung, U. Bunge, M. Schatz, and F. Thiele, *Restatement of the Spalart–Allmaras Eddy-Viscosity Model 
in Strain-Adaptive Formulation*, AIAA Journal, Vol. 41, no. 7, May 2012, [https://doi.org/10.2514/2.2089](https://arc.aiaa.org/doi/10.2514/2.2089)