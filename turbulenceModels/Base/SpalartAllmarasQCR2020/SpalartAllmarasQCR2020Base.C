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

#include "SpalartAllmarasQCR2020Base.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicSpalartAllmarasModel>
tmp<volScalarField>
SpalartAllmarasQCR2020Base<BasicSpalartAllmarasModel>::qcr2020Fw
(
    const volTensorField& gradU
) const
{
    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    const volScalarField Sbar(this->fv2(chi, fv1)*this->nuTilda_/sqr(this->kappa_*this->y_));
    const volScalarField OmegaS
    (
        sqrt(magSqr(skew(gradU)) + magSqr(symm(gradU)))
    );

    const dimensionedScalar epsS(OmegaS.dimensions(), SMALL);
    const volScalarField Slimit(Sbar + Cqcr2_*OmegaS);
    const volScalarField Stilda
    (
        pos(Slimit)*(OmegaS + Sbar)
      + neg(Slimit)
       *(
            OmegaS
          + OmegaS*(sqr(Cqcr2_)*OmegaS + Cqcr3_*Sbar)
           /max((Cqcr3_ - 2*Cqcr2_)*OmegaS - Sbar, epsS)
        )
    );

    const dimensionedScalar eps(Stilda.dimensions(), SMALL);
    volScalarField r
    (
        min
        (
            this->nuTilda_/(max(Stilda, eps)*sqr(this->kappa_*this->y_)),
            scalar(10)
        )
    );
    r.boundaryFieldRef() == 0;

    const volScalarField g(r + this->Cw2_*(pow6(r) - r));

    return g*pow((1 + pow6(this->Cw3_))/(pow6(g) + pow6(this->Cw3_)), 1.0/6.0);
}


template<class BasicSpalartAllmarasModel>
tmp<volSymmTensorField>
SpalartAllmarasQCR2020Base<BasicSpalartAllmarasModel>::qcr2020NonlinearStress
(
    const volVectorField& U
) const
{
    const volTensorField gradU(fvc::grad(U));
    const volSymmTensorField Rlin(-this->nut_*devTwoSymm(gradU));
    const dimensionedScalar epsGrad(gradU.dimensions(), SMALL);
    const volTensorField O(2*skew(gradU)/max(sqrt(magSqr(gradU)), epsGrad));
    const volScalarField fw(this->qcr2020Fw(gradU));
    const volScalarField Ccr1(Ccr1_*(scalar(1) + Cfw1_*fw));
    const volScalarField Ccr2(Ccr2_*(scalar(1) + Cfw2_*fw));

    return tmp<volSymmTensorField>::New
    (
        IOobject::groupName("qcr2020NonlinearStress", this->alphaRhoPhi_.group()),
      - Ccr1*twoSymm(O & Rlin) + Ccr2*this->nut_*this->Omega(gradU)*I
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicSpalartAllmarasModel>
SpalartAllmarasQCR2020Base<BasicSpalartAllmarasModel>::
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
)
:
    BasicSpalartAllmarasModel
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Ccr1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ccr1",
            this->coeffDict_,
            0.20
        )
    ),
    Ccr2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ccr2",
            this->coeffDict_,
            2.150537634408602
        )
    ),
    Cfw1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cfw1",
            this->coeffDict_,
            2.0
        )
    ),
    Cfw2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cfw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cqcr2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cqcr2",
            this->coeffDict_,
            0.7
        )
    ),
    Cqcr3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cqcr3",
            this->coeffDict_,
            0.9
        )
    )
{
    Info<< "QCR2020 stress correction: active" << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicSpalartAllmarasModel>
bool SpalartAllmarasQCR2020Base<BasicSpalartAllmarasModel>::read()
{
    if (BasicSpalartAllmarasModel::read())
    {
        Ccr1_.readIfPresent(this->coeffDict());
        Ccr2_.readIfPresent(this->coeffDict());
        Cfw1_.readIfPresent(this->coeffDict());
        Cfw2_.readIfPresent(this->coeffDict());
        Cqcr2_.readIfPresent(this->coeffDict());
        Cqcr3_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicSpalartAllmarasModel>
tmp<volSymmTensorField>
SpalartAllmarasQCR2020Base<BasicSpalartAllmarasModel>::R() const
{
    tmp<volSymmTensorField> tR(this->qcr2020NonlinearStress(this->U_));
    tR.ref() += -this->nut_*devTwoSymm(fvc::grad(this->U_));

    return tR;
}


template<class BasicSpalartAllmarasModel>
tmp<volSymmTensorField>
SpalartAllmarasQCR2020Base<BasicSpalartAllmarasModel>::devRhoReff() const
{
    return devRhoReff(this->U_);
}


template<class BasicSpalartAllmarasModel>
tmp<volSymmTensorField>
SpalartAllmarasQCR2020Base<BasicSpalartAllmarasModel>::devRhoReff
(
    const volVectorField& U
) const
{
    tmp<volSymmTensorField> tdevRhoReff
    (
        BasicSpalartAllmarasModel::devRhoReff(U)
    );

    tdevRhoReff.ref() += this->rho_*this->qcr2020NonlinearStress(U);

    return tdevRhoReff;
}


template<class BasicSpalartAllmarasModel>
tmp<fvVectorMatrix>
SpalartAllmarasQCR2020Base<BasicSpalartAllmarasModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
        fvc::div(this->rho_*this->qcr2020NonlinearStress(U))
      + BasicSpalartAllmarasModel::divDevRhoReff(U)
    );
}


template<class BasicSpalartAllmarasModel>
tmp<fvVectorMatrix>
SpalartAllmarasQCR2020Base<BasicSpalartAllmarasModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return
    (
        fvc::div(rho*this->qcr2020NonlinearStress(U))
      + BasicSpalartAllmarasModel::divDevRhoReff(rho, U)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
