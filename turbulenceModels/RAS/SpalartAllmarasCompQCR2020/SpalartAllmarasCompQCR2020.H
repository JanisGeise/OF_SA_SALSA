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

#ifndef Foam_SpalartAllmarasCompQCR2020_H
#define Foam_SpalartAllmarasCompQCR2020_H

#include "SpalartAllmarasComp.H"
#include "SpalartAllmarasQCR2020Base.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class SpalartAllmarasCompQCR2020
:
    public SpalartAllmarasQCR2020Base
    <
        SpalartAllmarasComp<BasicTurbulenceModel>
    >
{
    SpalartAllmarasCompQCR2020
    (
        const SpalartAllmarasCompQCR2020&
    ) = delete;

    void operator=(const SpalartAllmarasCompQCR2020&) = delete;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    TypeName("SpalartAllmarasCompQCR2020");

    SpalartAllmarasCompQCR2020
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

    virtual ~SpalartAllmarasCompQCR2020() = default;
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "SpalartAllmarasCompQCR2020.C"
#endif

#endif


// ************************************************************************* //
