# SALSA turbulence model for OpenFOAM
Implementation of the strain-adaptive linear Spalart-Allmaras (SALSA) turbulence model in OpenFOAM.

The repositories of [TUFRG](https://github.com/TUFRG/SAH-RANS-OF) and 
[mAlletto](https://gitlab.com/mAlletto/openfoamtutorials/-/tree/master/SpalartAllmarasRCsend)
were used as template for the structure of this implementation.

## Installation

To compile the SALSA turbulence model, execute `$./Allwmake`.

Once compiled, you should see `libSASALSAIncompressibleTurbulenceModel.so` `libSASALSACompressibleTurbulenceModel.so`
in `$FOAM_USER_LIBBIN`.

To remove the SALSA model, just execute the `$./Allwclean`.

## Setting up a case
To run a simulation with the SALSA model, first add the path to the library to the `controlDict`:

`libs			("libSASALSACompressibleTurbulenceModel.so");`

Then choose the correct turbulence model in the `turbulenceProperties`:

```
RAS
{
    RASModel            SpalartAllmarasSALSA;
    turbulence          on;
    printCoeffs         on;
}
```

There are three optional parameters `rhoInf, useRmod, useSmod`:
```
RAS
{
    RASModel            SpalartAllmarasSALSA;
    turbulence          on;
    printCoeffs         on;

    rhoInf              0.957837;   // required for useRmod = true
    useRmod             true;       // Edward's modification
    useSmod             true;      // strain rate instead of vorticity for Stilda
}
```
The parameter `rhoInf` denotes the free-stream densitiy $\rho_\infty$ and is only required for `useRmod   true`.
The parameter `useRmod` activates the Edward's modification, the parameter `useSmod` sets `STilda`to be the strain rate 
instead of the vorticity. Setting both paramerters to `true` results in the original SALSA turbulence model, linked in the references below.

You can find an example setup using the SALSA model [here](https://github.com/AndreWeiner/buffet_oat15).
## still TODO
- check if we can implement this more efficiently

## References
- template / directory structure taken from [TUFRG](https://github.com/TUFRG/SAH-RANS-OF) and 
[mAlletto](https://gitlab.com/mAlletto/openfoamtutorials/-/tree/master/SpalartAllmarasRCsend)
- T. Rung, U. Bunge, M. Schatz, and F. Thiele, *Restatement of the Spalart–Allmaras Eddy-Viscosity Model 
in Strain-Adaptive Formulation*, AIAA Journal, Vol. 41, no. 7, May 2012, [https://doi.org/10.2514/2.2089](https://arc.aiaa.org/doi/10.2514/2.2089)
- D.-M. Zimmermann, R. Mayer, T. Lutz, and E. Krämer, *Impact of model parameters of SALSA turbulence model
on transonic buffet prediction,* AIAA Journal, Vol. 56, no. 2, pp. 874–877, December 2017,
[https://arc.aiaa.org/doi/10.2514/1.J056193](https://arc.aiaa.org/doi/10.2514/1.J056193)