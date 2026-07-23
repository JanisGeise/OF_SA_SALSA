# SALSA turbulence models for OpenFOAM
Implementation of the strain-adaptive linear Spalart-Allmaras (SALSA) turbulence model in OpenFOAM.
The repository provides:

- `SpalartAllmarasSALSA`: the RANS model.
- `SpalartAllmarasSALSADDES`: a conservative DDES model using SALSA as the RANS baseline.
- `SpalartAllmarasQCR2020`: the standard SA RANS model with the QCR2020 nonlinear stress relation.
- `SpalartAllmarasComp`: standard SA with the mixing-layer compressibility correction.
- `SpalartAllmarasCompQCR2020`: SA-comp with the QCR2020 nonlinear stress relation.
- `SpalartAllmarasDDESQCR2020`: the standard SA-DDES model with the QCR2020 nonlinear stress relation.

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

To use the standard Spalart-Allmaras QCR2020 model, select
`SpalartAllmarasQCR2020`. This model is independent of SALSA and uses the
standard OpenFOAM SA transport equation with the QCR2020 nonlinear constitutive
stress relation.

```foam
RAS
{
    RASModel            SpalartAllmarasQCR2020;
    turbulence          on;
    printCoeffs         on;

    SpalartAllmarasQCR2020Coeffs
    {
        // QCR2020 defaults from Rumsey et al. (2020)
        Ccr1             0.20;
        Ccr2             2.150537634408602;
        Cfw1             2.0;
        Cfw2             0.3;
        Cqcr2            0.7;
        Cqcr3            0.9;

        // OpenFOAM's SA default is SA-noft2; set true for standard SA.
        ft2              false;
    }
}
```

For compressible mixing layers, `SpalartAllmarasComp` adds the SA-comp
destruction term of Spalart (2000) to the standard SA transport equation. It is
available only from the compressible library. The default coefficient is
`C5 3.5`:

```foam
RAS
{
    RASModel            SpalartAllmarasComp;
    turbulence          on;
    printCoeffs         on;

    SpalartAllmarasCompCoeffs
    {
        C5               3.5;
        ft2              false;
    }
}
```

To combine that transport correction with the QCR2020 stress relation, select
`SpalartAllmarasCompQCR2020`. Its coefficient dictionary accepts both `C5` and
the QCR2020 coefficients:

```foam
RAS
{
    RASModel            SpalartAllmarasCompQCR2020;
    turbulence          on;
    printCoeffs         on;

    SpalartAllmarasCompQCR2020Coeffs
    {
        C5               3.5;
        Ccr1             0.20;
        Ccr2             2.150537634408602;
        Cfw1             2.0;
        Cfw2             0.3;
        Cqcr2            0.7;
        Cqcr3            0.9;
        ft2              false;
    }
}
```

The implementation evaluates the local perfect-gas speed of sound as
$a^2=\gamma p/\rho$, consistently with the assumptions stated by NASA/TMR for
compressible SA. Setting `C5 0` removes the additional SA-comp term.

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

For the standard SA-DDES form with the QCR2020 stress relation, select
`SpalartAllmarasDDESQCR2020`. The same standard DDES and QCR2020 coefficients
can be supplied in `SpalartAllmarasDDESQCR2020Coeffs`.

### QCR2020 compatibility

The QCR2020 equations were taken from the NASA/TMR Spalart-Allmaras model page:
the stress relation, the wall-function-dependent `Ccr1''` and `Ccr2''`
coefficients, the normalized rotation tensor, and the smoother QCR2020
`Omega_s`-based `f_w` path are all directly extractable there. No conflicting
QCR2020 equations were found on that page.

QCR2020 is implemented here for standard SA RANS, SA-comp RANS, and standard
SA-DDES, but not as a SALSA option. SA-comp and QCR2020 are compatible: the
former modifies the transported-variable equation and the latter modifies the
constitutive stress relation. QCR2020 should not be combined with another QCR
variant. If a
future base model includes a true modeled turbulent kinetic energy term in the
Boussinesq relation, the QCR2020 isotropic approximation should be omitted to
avoid double counting.

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
- C. L. Rumsey, J.-R. Carlson, T. H. Pulliam, and P. R. Spalart,
*Improvements to the Quadratic Constitutive Relation Based on NASA Juncture Flow Data,*
AIAA Journal, Vol. 58, no. 10, pp. 4374-4384, 2020,
[https://doi.org/10.2514/1.J059683](https://doi.org/10.2514/1.J059683)
- P. R. Spalart, *Trends in Turbulence Treatments,* AIAA 2000-2306, June 2000,
[https://doi.org/10.2514/6.2000-2306](https://doi.org/10.2514/6.2000-2306)
