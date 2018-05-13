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

#include "opticallyThin.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "fvmSup.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"


using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(opticallyThin, 0);
        addToRadiationRunTimeSelectionTables(opticallyThin);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::opticallyThin::opticallyThin(const volScalarField& T)
:
    radiationModel(typeName, T),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    ),
    Tb_
    (
        dimensionedScalar("Tb", T.dimensions(), 300.0)
    )
{}


Foam::radiation::opticallyThin::opticallyThin(const dictionary& dict, const volScalarField& T)
:
    radiationModel(typeName, dict, T),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    ),
    Tb_
    (
        dimensionedScalar("Tb", T.dimensions(), 300.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::opticallyThin::~opticallyThin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::opticallyThin::read()
{
    if (radiationModel::read())
    {
        // nothing to read
        coeffs_.readIfPresent("Tb", Tb_); 
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::opticallyThin::calculate()
{}


Foam::tmp<Foam::volScalarField> Foam::radiation::opticallyThin::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            4.0*absorptionEmission_->eCont()*physicoChemical::sigma
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::opticallyThin::Ru() const
{
    const DimensionedField<scalar, volMesh> E =
        absorptionEmission_->ECont()().dimensionedInternalField();

    return  4.0*E;
}

Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::opticallyThin::Sh
(
    fluidThermo& thermo
) const
{
    volScalarField& he = thermo.he();
    const volScalarField T4(pow4(T_));
    
   
    return fvm::Su
    (
        -Rp()*( T4 - pow4(Tb_)),
        he        
    );
}
// ************************************************************************* //
