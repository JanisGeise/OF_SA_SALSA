/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "SpalartAllmarasSALSABase.H"
#include "wallDist.H"
#include "bound.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::chi() const
{
    return nuTilda_/this->nu();
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3("chi3", pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return scalar(1) - chi/(scalar(1) + chi*fv1);
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::ft2
(
    const volScalarField& chi
) const
{
    if (ft2_)
    {
        return Ct3_*exp(-Ct4_*sqr(chi));
    }

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "ft2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, Zero)
    );
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::Omega
(
    const volTensorField& gradU
) const
{
    return sqrt(2.0)*mag(skew(gradU));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::r
(
    const volScalarField& nur,
    const volScalarField& Stilda,
    const volScalarField& dTilda
) const
{
    const dimensionedScalar eps(Stilda.dimensions(), SMALL);

    tmp<volScalarField> tr =
        min(nur/(max(Stilda, eps)*sqr(kappa_*dTilda)), scalar(10));

    tr.ref().boundaryFieldRef() == 0;

    return tr;
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::rMod
(
    const volScalarField& nur,
    const volScalarField& Stilda,
    const volScalarField& dTilda
) const
{
    const dimensionedScalar eps(Stilda.dimensions(), SMALL);
    // const dimensionedScalar epsRho(this->rho_.dimensions(), SMALL);

    // compute
    // TODO: issue for rho = 0? shouldn't happen but who knows...
    tmp<volScalarField> Psi = sqrt(rhoInf_/this->rho_) * nuTilda_ / (sqr(kappa_) * sqr(dTilda));
    // doesn't compile: tmp<volScalarField> Psi = sqrt(rhoInf_/max(this->rho_, epsRho)) * nuTilda_ / (sqr(kappa_) * sqr(dTilda));

    // compute modified r, keep the limiter
    tmp<volScalarField> tr = min(1.6 * tanh(0.7 * Psi / max(Stilda, eps)), scalar(10));

    tr.ref().boundaryFieldRef() == 0;

    return tr;
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::fw
(
    const volScalarField& Stilda,
    const volScalarField& dTilda
) const
{
    // TODO: check how to make this more efficient
    if (rMod_)
    {
        const volScalarField::Internal r(this->rMod(nuTilda_, Stilda, dTilda)()());
        const volScalarField::Internal g(r + Cw2_*(pow6(r) - r));
        return g*pow((1 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
    }
    else
    {
        const volScalarField::Internal r(this->r(nuTilda_, Stilda, dTilda)()());
        const volScalarField::Internal g(r + Cw2_*(pow6(r) - r));
        return g*pow((1 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
    }
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU,
    const volScalarField& dTilda
) const
{
    const volScalarField Omega(this->Omega(gradU));

    return
        max
        (
            Omega + fv2(chi, fv1)*nuTilda_/sqr(kappa_*dTilda),
            Cs_*Omega
        );
}


template<class BasicEddyViscosityModel>
void SpalartAllmarasSALSABase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = nuTilda_*fv1;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicEddyViscosityModel>
void SpalartAllmarasSALSABase<BasicEddyViscosityModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}

template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::GammaEff
(
    const volScalarField& Stilda,
    const volScalarField& dTilda
) const
{
    const volScalarField::Internal r(this->r(nuTilda_, Stilda, dTilda)()());

    const volScalarField::Internal alpha1
    (
        pow(scalar(1.01)*r, scalar(0.65))
    );

    const volScalarField::Internal alpha2
    (
        pow
        (
            max
            (
                scalar(0),
                scalar(1) - tanh(this->chi()/scalar(68))
            ),
            scalar(0.65)
        )
    );

    const volScalarField::Internal gamma(max(alpha1, alpha2));

    return
        sqrt(
            min(
                scalar(1.25),
                max(gamma, scalar(0.75))
            )
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
SpalartAllmarasSALSABase<BasicEddyViscosityModel>::SpalartAllmarasSALSABase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
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

    sigmaNut_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaNut",
            this->coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    Cb1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb1",
            this->coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb2",
            this->coeffDict_,
            0.0
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),        // TODO: leave that in here even though not used anymore?
    Cw2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cv1",
            this->coeffDict_,
            7.1
        )
    ),
    Cs_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.3
        )
    ),
    ck_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "ck",
            this->coeffDict_,
            0.07
        )
    ),
    ft2_
    (
        Switch::getOrAddToDict
        (
            "ft2",
            this->coeffDict_,
            false
        )
    ),
    Ct3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ct3",
            this->coeffDict_,
            1.2
        )
    ),
    Ct4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ct4",
            this->coeffDict_,
            0.5
        )
    ),

    rhoInf_
    (
        dimensionedScalar
        (
            "rhoInf",
            dimDensity,
            dimensioned<scalar>::getOrAddToDict
            (
                "rhoInf",
                this->coeffDict_,
                1.0
            ).value()
        )
    ),

    rMod_
    (
        Switch::getOrAddToDict
        (
            "useRmod",
            this->coeffDict_,
            false
        )
    ),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    y_(wallDist::New(this->mesh_).y())
{
    if (rMod_)
    {
        Info<< "modified r term: active" << nl;
    }
    else
    {
        Info<< "modified r term: inactive" << nl;
    }

    if (ft2_)
    {
        Info<< "ft2 term: active" << nl;
    }
    else
    {
        Info<< "ft2 term: inactive" << nl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
bool SpalartAllmarasSALSABase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;       // TODO: leave that in here even though not used anymore?
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        ck_.readIfPresent(this->coeffDict());

        ft2_.readIfPresent("ft2", this->coeffDict());
        Ct3_.readIfPresent(this->coeffDict());
        Ct4_.readIfPresent(this->coeffDict());
        rhoInf_.readIfPresent(this->coeffDict());

        if (rMod_)
        {
            Info<< "    modified r term: active" << nl;
        }
        else
        {
            Info<< "    modified r term: inactive" << nl;
        }

        if (ft2_)
        {
            Info<< "    ft2 term: active" << nl;
        }
        else
        {
            Info<< "    ft2 term: inactive" << nl;
        }

        return true;
    }

    return false;
}


template<class BasicEddyViscosityModel>
tmp<volScalarField>
SpalartAllmarasSALSABase<BasicEddyViscosityModel>::DnuTildaEff() const
{
    return tmp<volScalarField>::New
    (
        IOobject::groupName("DnuTildaEff", this->alphaRhoPhi_.group()),
        (nuTilda_ + this->nu())/sigmaNut_
    );
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::k() const
{
    // (B:Eq. 4.50)
    const scalar Cmu = 0.09;
    const auto fv1 = this->fv1(chi());

    return tmp<volScalarField>::New
    (
        IOobject::groupName("k", this->alphaRhoPhi_.group()),
        cbrt(fv1)*nuTilda_*::sqrt(scalar(2)/Cmu)*mag(symm(fvc::grad(this->U_)))
    );
}

template<class BasicEddyViscosityModel>
tmp<volScalarField>
SpalartAllmarasSALSABase<BasicEddyViscosityModel>::epsilon() const
{
    // (B:Eq. 4.50)
    const scalar Cmu = 0.09;
    const auto fv1 = this->fv1(chi());
    const dimensionedScalar nutSMALL(nuTilda_.dimensions(), SMALL);

    return tmp<volScalarField>::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        sqrt(fv1)*sqr(::sqrt(Cmu)*this->k())/(nuTilda_ + this->nut_ + nutSMALL)
    );
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> SpalartAllmarasSALSABase<BasicEddyViscosityModel>::omega() const
{
    // (P:p. 384)
    const scalar betaStar = 0.09;
    const dimensionedScalar k0(sqr(dimLength/dimTime), SMALL);

    return tmp<volScalarField>::New
    (
        IOobject::groupName("omega", this->alphaRhoPhi_.group()),
        this->epsilon()/(betaStar*(this->k() + k0))
    );
}


template<class BasicEddyViscosityModel>
void SpalartAllmarasSALSABase<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    {
        // Local references
        const alphaField& alpha = this->alpha_;
        const rhoField& rho = this->rho_;
        const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
        const volVectorField& U = this->U_;
        fv::options& fvOptions(fv::options::New(this->mesh_));

        BasicEddyViscosityModel::correct();

        const volScalarField chi(this->chi());
        const volScalarField fv1(this->fv1(chi));
        const volScalarField ft2(this->ft2(chi));

        tmp<volTensorField> tgradU = fvc::grad(U);
        volScalarField dTilda(this->dTilda(chi, fv1, tgradU()));
        volScalarField Stilda(this->Stilda(chi, fv1, tgradU(), dTilda));
        tgradU.clear();

        // compute Cb1*sqrt(gamma) for SALSA modification and update Cw1 accordingly
        const volScalarField::Internal Cb1Eff = Cb1_ * GammaEff(Stilda, dTilda)();
        const volScalarField::Internal Cw1Eff = Cb1Eff/sqr(kappa_) + (scalar(1) + Cb2_)/sigmaNut_;

        tmp<fvScalarMatrix> nuTildaEqn
        (
            fvm::ddt(alpha, rho, nuTilda_)
          + fvm::div(alphaRhoPhi, nuTilda_)
          - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
          - Cb2_/sigmaNut_*alpha()*rho()*magSqr(fvc::grad(nuTilda_)()())
         ==
            Cb1Eff*alpha()*rho()*Stilda()*nuTilda_()*(scalar(1) - ft2())
          - fvm::Sp
            (
                (Cw1Eff*fw(Stilda, dTilda) - Cb1Eff/sqr(kappa_)*ft2())
               *alpha()*rho()*nuTilda_()/sqr(dTilda()),
                nuTilda_
            )
          + fvOptions(alpha, rho, nuTilda_)
        );

        nuTildaEqn.ref().relax();
        fvOptions.constrain(nuTildaEqn.ref());
        solve(nuTildaEqn);
        fvOptions.correct(nuTilda_);
        bound(nuTilda_, dimensionedScalar(nuTilda_.dimensions(), Zero));
        nuTilda_.correctBoundaryConditions();
    }

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
