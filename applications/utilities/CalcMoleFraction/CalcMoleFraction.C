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

Application
    CalcTransport

Description
    Utility for writing out the mole fraction fields from mass fractions.

\*---------------------------------------------------------------------------*/
#include "calc.H"
#include "fvc.H"
#include "fvCFD.H"
#include "turbulenceModel.H"
#include "psiCombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "moleFraction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    autoPtr<combustionModels::psiCombustionModel> reaction
    (
        combustionModels::psiCombustionModel::New(mesh)
    );

    psiReactionThermo& thermo = reaction->thermo();
    thermo.validate(args.executable(), "h", "e");

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();
    moleFraction moleFraction_(composition, mesh);
    PtrList<volScalarField>& X = moleFraction_.X();

    forAll(X, i)
    {
        X[i].write();
    }
    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
