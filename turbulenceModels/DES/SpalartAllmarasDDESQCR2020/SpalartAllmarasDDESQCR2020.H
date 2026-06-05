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

#ifndef Foam_SpalartAllmarasDDESQCR2020_H
#define Foam_SpalartAllmarasDDESQCR2020_H

#include "SpalartAllmarasDES.H"
#include "SpalartAllmarasDDES.H"
#include "SpalartAllmarasQCR2020Base.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                   Class SpalartAllmarasDDESQCR2020 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasDDESQCR2020
:
    public SpalartAllmarasQCR2020Base
    <
        SpalartAllmarasDDES<BasicTurbulenceModel>
    >
{
    // Private Member Functions

        //- No copy construct
        SpalartAllmarasDDESQCR2020
        (
            const SpalartAllmarasDDESQCR2020&
        ) = delete;

        //- No copy assignment
        void operator=(const SpalartAllmarasDDESQCR2020&) = delete;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmarasDDESQCR2020");


    // Constructors

        //- Construct from components
        SpalartAllmarasDDESQCR2020
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


    //- Destructor
    virtual ~SpalartAllmarasDDESQCR2020() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasDDESQCR2020.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
