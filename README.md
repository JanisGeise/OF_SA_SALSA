# SALSA turbulence models for OpenFOAM
Implementation of the strain-adaptive linear Spalart-Allmaras (SALSA) turbulence model in OpenFOAM.
The repository provides:

- `SpalartAllmarasSALSA`: the RANS model.
- `SpalartAllmarasSALSADDES`: a conservative DDES model using SALSA as the RANS baseline.

The repositories of [TUFRG](https://github.com/TUFRG/SAH-RANS-OF) and 
[mAlletto](https://gitlab.com/mAlletto/openfoamtutorials/-/tree/master/SpalartAllmarasRCsend)
were used as template for the structure of this implementation.

## Installation

To compile the SALSA turbulence models, execute:

```bash
./Allwmake
```

Once compiled, you should see `libSASALSAIncompressibleTurbulenceModel.so` and
`libSASALSACompressibleTurbulenceModel.so` in `$FOAM_USER_LIBBIN`.

To remove the compiled model libraries and build files, execute:

```bash
./Allwclean
```

## Setting up a case
To run a simulation with one of the SALSA models, first add the corresponding library to the `controlDict`.
For compressible solvers:

```foam
libs    ("libSASALSACompressibleTurbulenceModel.so");
```

For incompressible solvers:

```foam
libs    ("libSASALSAIncompressibleTurbulenceModel.so");
```

### RANS model

Choose the RANS model in `constant/turbulenceProperties`:

```foam
simulationType  RAS;

RAS
{
    RASModel            SpalartAllmarasSALSA;
    turbulence          on;
    printCoeffs         on;
}
```

There are three optional parameters: `rhoInf`, `useRmod`, and `useSmod`.

```foam
RAS
{
    RASModel            SpalartAllmarasSALSA;
    turbulence          on;
    printCoeffs         on;

    rhoInf              0.957837;   // required for useRmod = true
    useRmod             true;       // Edward's modification
    useSmod             true;       // strain rate instead of vorticity for Stilda
}
```
The parameter `rhoInf` denotes the free-stream density $\rho_\infty$ and is only required for `useRmod true`.
The parameter `useRmod` activates Edwards' modification, and `useSmod` sets `Stilda` to be the strain rate
instead of the vorticity. Setting both parameters to `true` results in the original SALSA turbulence model,
linked in the references below.

### DDES model

The `SpalartAllmarasSALSADDES` model is selected as an LES model:

```foam
simulationType  LES;

LES
{
    LESModel            SpalartAllmarasSALSADDES;
    turbulence          on;
    printCoeffs         on;

    delta               cubeRootVol;
}
```

The model uses the SALSA transport equation and source/destruction terms, but replaces the RANS wall-distance
length scale with the DDES length scale. Optional DDES and SALSA coefficients can be supplied in the
`SpalartAllmarasSALSADDESCoeffs` dictionary:

```foam
LES
{
    LESModel            SpalartAllmarasSALSADDES;
    turbulence          on;
    printCoeffs         on;

    delta               cubeRootVol;

    SpalartAllmarasSALSADDESCoeffs
    {
        // DDES settings
        CDES              0.65;
        shielding         standard;     // standard or ZDES2020
        lowReCorrection   true;

        // Only used with shielding ZDES2020
        Cd3               25;
        Cd4               0.03;
        betaZDES          2.5;
        usefP2            false;

        // SALSA settings
        rhoInf            0.957837;     // required for useRmod = true
        useRmod           true;
        useSmod           true;
    }
}
```

Note on `useSigma`: OpenFOAM's standard `SpalartAllmarasDDES` has a `useSigma` switch for a sigma-based
grey-area enhancement. This implementation deliberately does not expose `useSigma`, because it changes the
production measure and can conflict conceptually with SALSA's own `useSmod` strain-rate formulation. The first
`SpalartAllmarasSALSADDES` version is therefore conservative: SALSA source/destruction terms plus DDES `dTilda`
and standard or `ZDES2020` shielding.

The source-code consistent equations for the standard SA, SALSA, and DDES variants are collected in
[equations.md](equations.md).

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
