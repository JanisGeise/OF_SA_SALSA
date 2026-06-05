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

#ifndef Foam_SpalartAllmarasQCR2020Base_H
#define Foam_SpalartAllmarasQCR2020Base_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class SpalartAllmarasQCR2020Base Declaration
\*---------------------------------------------------------------------------*/

template<class BasicSpalartAllmarasModel>
class SpalartAllmarasQCR2020Base
:
    public BasicSpalartAllmarasModel
{
    // Private Member Functions

        //- No copy construct
        SpalartAllmarasQCR2020Base
        (
            const SpalartAllmarasQCR2020Base&
        ) = delete;

        //- No copy assignment
        void operator=(const SpalartAllmarasQCR2020Base&) = delete;


protected:

    // Protected Data

        // Model constants

            dimensionedScalar Ccr1_;
            dimensionedScalar Ccr2_;
            dimensionedScalar Cfw1_;
            dimensionedScalar Cfw2_;
            dimensionedScalar Cqcr2_;
            dimensionedScalar Cqcr3_;


    // Protected Member Functions

        tmp<volScalarField> qcr2020Fw(const volTensorField& gradU) const;

        tmp<volSymmTensorField> qcr2020NonlinearStress
        (
            const volVectorField& U
        ) const;


public:

    typedef typename BasicSpalartAllmarasModel::alphaField alphaField;
    typedef typename BasicSpalartAllmarasModel::rhoField rhoField;
    typedef typename BasicSpalartAllmarasModel::transportModel transportModel;


    // Constructors

        //- Construct from components
        SpalartAllmarasQCR2020Base
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName,
            const word& type
        );


    //- Destructor
    virtual ~SpalartAllmarasQCR2020Base() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the effective stress tensor based on a given velocity field
        virtual tmp<volSymmTensorField> devRhoReff
        (
            const volVectorField& U
        ) const;

        //- Return the divergence of the effective stress tensor
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Return the divergence of the effective stress tensor
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasQCR2020Base.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
