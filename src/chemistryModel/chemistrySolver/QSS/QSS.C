/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | Code based on OpenFOAM
  \\    /   O peration     |
   \\  /    A nd           | Copyright (C) Adhiraj Dasgupta
    \\/     M anipulation  |                     
-------------------------------------------------------------------------------
 License
     This file is a derivative work of OpenFOAM.
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

#include "QSS.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::QSS<ChemistryModel>::QSS
(
    const fvMesh& mesh
)
:
    chemistrySolver<ChemistryModel>(mesh),
    coeffsDict_(this->subDict("QSSCoeffs")),
    odeSolver_(ODESolver::New(*this, coeffsDict_)),
    cTp_(this->nEqns())
{
    Ytmp0 = new double [this->nSpecie()];
    dYdt = new double [this->nSpecie()];
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::QSS<ChemistryModel>::~QSS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::QSS<ChemistryModel>::solve
(
    scalarField& c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{

    label nSpecie = this->nSpecie();

    // Copy the concentration, T and P to the total solve-vector
    for (register int i=0; i<nSpecie; i++)
    {
        cTp_[i] = c[i];
    }
    cTp_[nSpecie] = T;
    cTp_[nSpecie+1] = p;

    odeSolver_->solve(0, deltaT, cTp_, subDeltaT);

    for (register int i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, cTp_[i]);
    }
    T = cTp_[nSpecie];
    p = cTp_[nSpecie+1];
}


template<class ChemistryModel>
void Foam::QSS<ChemistryModel>::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{

    label nSpecie = this->nSpecie();
    double T = c[nSpecie];
    const scalar p = c[nSpecie + 1];
    //Convert Pa to dynes/cm2
    double pa = p*10.0;

    // constant pressure
    // dT/dt = ...
    scalar rho = 0.0;
    scalar cSum = 0.0;
    for (label i = 0; i < nSpecie; i++)
    {
        const scalar W = this->specieThermo_[i].W();
        cSum += c[i];
        rho += W*c[i];
    }
    
    // Calculate the mass fractions since the routine needs Y
    for (label i = 0; i < nSpecie; i++)
    {
        const scalar W = this->specieThermo_[i].W();
        Ytmp0[i] = c[i]*W/rho;        
    }    
    ckwyp_(&pa, &T, Ytmp0, NULL, NULL, dYdt);
    for (label i = 0; i < nSpecie; i++)
    {
        dcdt[i] = dYdt[i]*1e3;
    }
        
    scalar cp = 0.0;
    for (label i=0; i<nSpecie; i++)
    {
        cp += c[i]*(this->specieThermo_[i].cp(p, T));
    }
    cp /= rho;

    scalar dT = 0.0;
    for (label i = 0; i < nSpecie; i++)
    {
        const scalar hi = this->specieThermo_[i].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    dcdt[nSpecie] = -dT;

    // dp/dt = ...
    dcdt[nSpecie + 1] = 0.0;    
}

template<class ChemistryModel>
void Foam::QSS<ChemistryModel>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    //- This uses a difference quotient Jacobian. The original C++
    //  code snippet was provided by Dr. Abdurrehman Imren in private
    //  communication.
    label nSpecie = this->nSpecie();
    label nEqns = this->nEqns();
 
    const scalar uround = 2.220446049250313E-016;
    scalarField c3 = c;
    scalarField dcdtH(nEqns, 0.0);
    derivatives(t, c, dcdt);
    for (label i = 0; i < nEqns; i++)
    {
        scalar delt = sqrt(uround*max(1.0e-5,mag(c[i])));
        c3[i] = c[i] + delt;
        derivatives(t, c3, dcdtH);
        for (label j = 0; j < nEqns; j++)
        {   
            dfdc[j][i] = (dcdtH[j] - dcdt[j])/delt;
        }
        c3[i] = c[i];
    }
    
}

template<class ChemistryModel>
Foam::tmp<Foam::scalarField> Foam::QSS<ChemistryModel>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    notImplemented
    (
        "QSS::omega"
        "("
            "const scalarField&, "
            "const scalar, "
            "const scalar"
        ") const"
    );  
}

template<class ChemistryModel>
template<class ThermoType>
Foam::scalar Foam::QSS<ChemistryModel>::omega
(
    const Reaction<ThermoType>& r,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef    
) const
{
    notImplemented
    (
        "QSS::omega"
        "("
             "const Reaction<ThermoType>&, "
             "const scalarField&, "
             "const scalar, "
             "const scalar, "
             "scalar&, "
             "scalar&, "
             "label&, "
             "scalar&, "
             "scalar&, "
             "label&"
        ") const"
    );
}

template<class ChemistryModel>
Foam::scalar Foam::QSS<ChemistryModel>::omegaI
(
    label iReaction,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef    
) const
{
    notImplemented
    (
        "QSS::omegaI"
        "("
             "label, "
             "const scalarField&, "
             "const scalar, "
             "const scalar, "
             "scalar&, "
             "scalar&, "
             "label&, "
             "scalar&, "
             "scalar&, "
             "label&"
        ") const"
    );
}

template<class ChemistryModel>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::QSS<ChemistryModel>::calculateRR
(
    const label reactionI,
    const label specieI    
) const
{
    notImplemented
    (
        "QSS::calculateRR"
        "("
             "const label, "
             "const label"
        ") const"    
    );
}

template<class ChemistryModel>
Foam::tmp<Foam::volScalarField> Foam::QSS<ChemistryModel>::tc() const
{
    notImplemented
    (
        "QSS::tc() const"
    );
}

template<class ChemistryModel>
void Foam::QSS<ChemistryModel>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }
    
    label nSpecie = this->nSpecie();
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        double Ti = T[celli];
        const scalar pi = p[celli];
        //- Convert Pa to dynes/cm2
        scalar pa = pi*10.0;

        scalarField c(nSpecie, 0.0);
        for (label i=0; i<nSpecie; i++)
        {
            Ytmp0[i] = this->Y_[i][celli];           
        }

        ckwyp_(&pa, &Ti, Ytmp0, NULL, NULL, dYdt);
        scalarField dcdt(nSpecie);
        for (label i = 0; i < nSpecie; i++)
        {
            dcdt[i] = dYdt[i]*1e3;
        }
        for (label i=0; i<nSpecie; i++)
        {
            this->RR_[i][celli] = dcdt[i]*this->specieThermo_[i].W();
        }
    }    
}
// ************************************************************************* //
