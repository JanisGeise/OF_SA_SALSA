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

\*---------------------------------------------------------------------------*/

#include "SpalartAllmarasSALSADDES.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
const Foam::Enum
<
    typename Foam::LESModels::SpalartAllmarasSALSADDES<BasicTurbulenceModel>::shieldingMode
>
Foam::LESModels::SpalartAllmarasSALSADDES<BasicTurbulenceModel>::shieldingModeNames
({
    { shieldingMode::standard, "standard" },
    { shieldingMode::ZDES2020, "ZDES2020" },
});


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasSALSADDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    const volScalarField r(this->r(this->nuEff(), magGradU, this->y_));

    tmp<volScalarField> tfd = 1 - tanh(pow(Cd1_*r, Cd2_));

    switch (shielding_)
    {
        case shieldingMode::standard:
        {
            return tfd;
        }
        case shieldingMode::ZDES2020:
        {
            auto maxEps = [](const volScalarField& fld, const scalar eps)
            {
                return max(fld, dimensionedScalar(fld.dimensions(), eps));
            };

            volScalarField& fdStd = tfd.ref();
            const auto& nuTilda = this->nuTilda_;
            const volVectorField& n = wallDist::New(this->mesh_).n();

            const volScalarField GnuTilda
            (
                Cd3_*maxEps(fvc::grad(nuTilda) & n, Zero)
              / (maxEps(magGradU, SMALL)*this->kappa_*this->y_)
            );

            volScalarField fdGnuTilda(1 - tanh(pow(Cd1_*GnuTilda, Cd2_)));
            const volScalarField GOmega
            (
              - (fvc::grad(mag(fvc::curl(this->U_))) & n)
              * sqrt(nuTilda/maxEps(pow3(magGradU), SMALL))
            );
            const volScalarField alpha((7.0/6.0*Cd4_ - GOmega)/(Cd4_/6.0));
            const volScalarField fRGOmega
            (
                pos(Cd4_ - GOmega)
              + 1.0
               /(1 + exp(min(-6*alpha/max(1 - sqr(alpha), SMALL), scalar(50))))
               *pos(4*Cd4_/3.0 - GOmega)*pos(GOmega - Cd4_)
            );

            if (usefP2_)
            {
                fdGnuTilda *=
                    (1.0 - tanh(pow(Cd1_*betaZDES_*r, Cd2_)))
                  / maxEps(fdStd, SMALL);
            }

            fdStd *= 1 - (1 - fdGnuTilda)*fRGOmega;
        }
    }

    return tfd;
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasSALSADDES<BasicTurbulenceModel>::psi
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    auto tpsi = tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::scopedName(this->type(), "psi"),
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    );

    if (lowReCorrection_)
    {
        auto& psi = tpsi.ref();

        const volScalarField fv2(this->fv2(chi, fv1));
        const volScalarField ft2(this->ft2(chi));

        psi =
            sqrt
            (
                min
                (
                    scalar(100),
                    (1 - this->Cb1_/(this->Cw1_*sqr(this->kappa_)*fwStar_)
                   *(ft2 + (1 - ft2)*fv2))
                   /max(SMALL, (fv1*max(scalar(1e-10), 1 - ft2)))
                )
            );
    }

    return tpsi;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasSALSADDES<BasicTurbulenceModel>::lengthScaleLES
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return psi(chi, fv1)*CDES_*this->delta();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasSALSADDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    const volScalarField& lRAS(this->y_);
    const volScalarField lLES(this->lengthScaleLES(chi, fv1));
    const dimensionedScalar l0(dimLength, Zero);

    return max
    (
        lRAS - fd(mag(gradU))*max(lRAS - lLES, l0),
        dimensionedScalar("small", dimLength, SMALL)
    );
}


template<class BasicTurbulenceModel>
void SpalartAllmarasSALSADDES<BasicTurbulenceModel>::correctNut()
{
    // Correct the turbulence viscosity
    SpalartAllmarasSALSABase<DESModel<BasicTurbulenceModel>>::correctNut();

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmarasSALSADDES<BasicTurbulenceModel>::SpalartAllmarasSALSADDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    SpalartAllmarasSALSABase<DESModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    shielding_
    (
        shieldingModeNames.getOrDefault
        (
            "shielding",
            this->coeffDict_,
            shieldingMode::standard
        )
    ),
    CDES_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CDES",
            this->coeffDict_,
            0.65
        )
    ),
    lowReCorrection_
    (
        Switch::getOrAddToDict
        (
            "lowReCorrection",
            this->coeffDict_,
            true
        )
    ),
    fwStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "fwStar",
            this->coeffDict_,
            0.424
        )
    ),
    Cd1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cd1",
            this->coeffDict_,
            8
        )
    ),
    Cd2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cd2",
            this->coeffDict_,
            3
        )
    ),
    Cd3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cd3",
            this->coeffDict_,
            25
        )
    ),
    Cd4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cd4",
            this->coeffDict_,
            0.03
        )
    ),
    betaZDES_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaZDES",
            this->coeffDict_,
            2.5
        )
    ),
    usefP2_
    (
        Switch::getOrAddToDict
        (
            "usefP2",
            this->coeffDict_,
            false
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);

        switch (shielding_)
        {
            case shieldingMode::standard:
            {
                Info<< "shielding function: standard DDES "
                    <<  "(Spalart et al., 2006)"
                    << nl;
                break;
            }
            case shieldingMode::ZDES2020:
            {
                Info<< "shielding function: ZDES mode 2 (Deck & Renard, 2020)"
                    << nl;
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unrecognised 'shielding' option: "
                    << shieldingModeNames[shielding_]
                    << exit(FatalError);
            }
        }

        if (usefP2_)
        {
            Info<< "fP2 term: active" << nl;
        }
        else
        {
            Info<< "fP2 term: inactive" << nl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasSALSADDES<BasicTurbulenceModel>::read()
{
    if (SpalartAllmarasSALSABase<DESModel<BasicTurbulenceModel>>::read())
    {
        shieldingModeNames.readIfPresent
        (
            "shielding",
            this->coeffDict(),
            shielding_
        );

        CDES_.readIfPresent(this->coeffDict());
        lowReCorrection_.readIfPresent("lowReCorrection", this->coeffDict());
        fwStar_.readIfPresent(this->coeffDict());
        Cd1_.readIfPresent(this->coeffDict());
        Cd2_.readIfPresent(this->coeffDict());
        Cd3_.readIfPresent(this->coeffDict());
        Cd4_.readIfPresent(this->coeffDict());
        betaZDES_.readIfPresent(this->coeffDict());
        usefP2_.readIfPresent("usefP2", this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasSALSADDES<BasicTurbulenceModel>::LESRegion() const
{
    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::scopedName("DDES", "LESRegion"),
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        neg(dTilda(chi, fv1, fvc::grad(this->U_)) - this->y_)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasSALSADDES<BasicTurbulenceModel>::fd() const
{
    return fd(mag(fvc::grad(this->U_)));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
