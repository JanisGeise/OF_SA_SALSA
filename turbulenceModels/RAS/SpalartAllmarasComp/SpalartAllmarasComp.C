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

#include "SpalartAllmarasComp.H"
#include "fvOptions.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
SpalartAllmarasComp<BasicTurbulenceModel>::SpalartAllmarasComp
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
    SpalartAllmaras<BasicTurbulenceModel>
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
    C5_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C5",
            this->coeffDict_,
            3.5
        )
    )
{
    Info<< "SA mixing-layer compressibility correction: active" << nl;

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


template<class BasicTurbulenceModel>
bool SpalartAllmarasComp<BasicTurbulenceModel>::read()
{
    if (SpalartAllmaras<BasicTurbulenceModel>::read())
    {
        C5_.readIfPresent(this->coeffDict());
        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void SpalartAllmarasComp<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    {
        const alphaField& alpha = this->alpha_;
        const rhoField& rho = this->rho_;
        const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
        const volVectorField& U = this->U_;
        fv::options& fvOptions(fv::options::New(this->mesh_));

        typedef eddyViscosity<RASModel<BasicTurbulenceModel>> evmType;
        evmType::correct();

        const volScalarField chi(this->chi());
        const volScalarField fv1(this->fv1(chi));
        const volScalarField ft2(this->ft2(chi));

        tmp<volTensorField> tgradU = fvc::grad(U);
        volScalarField dTilda(this->dTilda(chi, fv1, tgradU()));
        volScalarField Stilda(this->Stilda(chi, fv1, tgradU(), dTilda));

        // TMR SA-comp assumes a perfect gas: a^2 = gamma*p/rho.
        const volScalarField aSqr
        (
            this->transport().gamma()
           *this->transport().p()/rho
        );

        tmp<fvScalarMatrix> nuTildaEqn
        (
            fvm::ddt(alpha, rho, this->nuTilda_)
          + fvm::div(alphaRhoPhi, this->nuTilda_)
          - fvm::laplacian
            (
                alpha*rho*this->DnuTildaEff(),
                this->nuTilda_
            )
          - this->Cb2_/this->sigmaNut_*alpha()*rho()
           *magSqr(fvc::grad(this->nuTilda_)()())
         ==
            this->Cb1_*alpha()*rho()*Stilda()*this->nuTilda_()
           *(scalar(1) - ft2())
          - fvm::Sp
            (
                (this->Cw1_*this->fw(Stilda, dTilda)
               - this->Cb1_/sqr(this->kappa_)*ft2())
               *alpha()*rho()*this->nuTilda_()/sqr(dTilda()),
                this->nuTilda_
            )
          - fvm::Sp
            (
                C5_*alpha()*rho()*this->nuTilda_()
               *magSqr(tgradU()())/aSqr(),
                this->nuTilda_
            )
          + fvOptions(alpha, rho, this->nuTilda_)
        );

        tgradU.clear();

        nuTildaEqn.ref().relax();
        fvOptions.constrain(nuTildaEqn.ref());
        solve(nuTildaEqn);
        fvOptions.correct(this->nuTilda_);
        bound
        (
            this->nuTilda_,
            dimensionedScalar(this->nuTilda_.dimensions(), Zero)
        );
        this->nuTilda_.correctBoundaryConditions();
    }

    this->correctNut();
}

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
