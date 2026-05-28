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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::LESModels::SpalartAllmarasSALSADDES

Description
    Delayed detached-eddy simulation model using SALSA as the RANS baseline.

    The transport equation and SALSA source/destruction modifications are
    inherited from SpalartAllmarasSALSABase. The DES length-scale replacement
    and shielding functions follow the OpenFOAM v2506 SpalartAllmarasDDES
    implementation.

SourceFiles
    SpalartAllmarasSALSADDES.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SpalartAllmarasSALSADDES_H
#define Foam_SpalartAllmarasSALSADDES_H

#include "SpalartAllmarasSALSABase.H"
#include "LESeddyViscosity.H"
#include "DESModel.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class SpalartAllmarasSALSADDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasSALSADDES
:
    public SpalartAllmarasSALSABase<DESModel<BasicTurbulenceModel>>
{
public:

    // Public enumerations

        //- Shielding modes
        enum class shieldingMode
        {
            standard,
            ZDES2020
        };

        //- Shielding mode names
        static const Enum<shieldingMode> shieldingModeNames;


private:

    // Private Member Functions

        //- Return the shielding function
        tmp<volScalarField> fd(const volScalarField& magGradU) const;

        //- No copy construct
        SpalartAllmarasSALSADDES(const SpalartAllmarasSALSADDES&) = delete;

        //- No copy assignment
        void operator=(const SpalartAllmarasSALSADDES&) = delete;


protected:

    // Protected Data

        //- Shielding mode
        shieldingMode shielding_;

        // Model coefficients

            //- DES coefficient
            dimensionedScalar CDES_;

            //- Flag for low Reynolds number correction
            Switch lowReCorrection_;
            dimensionedScalar fwStar_;

            dimensionedScalar Cd1_;
            dimensionedScalar Cd2_;
            dimensionedScalar Cd3_;
            dimensionedScalar Cd4_;
            dimensionedScalar betaZDES_;
            Switch usefP2_;


    // Protected Member Functions

        //- Return the low Reynolds number correction function
        virtual tmp<volScalarField> psi
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        //- Return the LES length scale
        virtual tmp<volScalarField> lengthScaleLES
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        //- Return the length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU
        ) const;

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmarasSALSADDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasSALSADDES
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
    virtual ~SpalartAllmarasSALSADDES() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;

        //- Return the shielding function
        virtual tmp<volScalarField> fd() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasSALSADDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
