/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#ifndef Foam_SpalartAllmarasComp_H
#define Foam_SpalartAllmarasComp_H

#include "SpalartAllmaras.H"

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                    Class SpalartAllmarasComp Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasComp
:
    public SpalartAllmaras<BasicTurbulenceModel>
{
    // Private Member Functions

        SpalartAllmarasComp(const SpalartAllmarasComp&) = delete;
        void operator=(const SpalartAllmarasComp&) = delete;


protected:

    // Model constants

        //- Mixing-layer compressibility correction coefficient
        dimensionedScalar C5_;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    TypeName("SpalartAllmarasComp");

    SpalartAllmarasComp
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );

    virtual ~SpalartAllmarasComp() = default;

    virtual bool read();

    //- Correct nuTilda including the SA-comp destruction term
    virtual void correct();
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "SpalartAllmarasComp.C"
#endif

#endif

// ************************************************************************* //
