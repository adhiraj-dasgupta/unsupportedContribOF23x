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
    CalcMixtureFraction

Description
    Post-processing utility to calculate the mixture fraction based on Bilger's
    formula. Requires the mixtureFractionDict
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "dictionary.H"
#include "psiCombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "OSspecific.H"
#include "chemkinReader.H"
#include "makeGraph.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    
//  Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "setRootCase.H"
#   include "createTime.H"
    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"    
#   include "createFields.H"
    speciesTable species;
    IOdictionary mixtureFractionDict
    (
        IOobject
        (
            "mixtureFractionDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    dictionary fuel(mixtureFractionDict.subDict("fuel"));
    dictionary oxidizer(mixtureFractionDict.subDict("oxidizer"));
    const word fractionBasis(mixtureFractionDict.lookup("fractionBasis"));
    if ((fractionBasis != "mass") && (fractionBasis != "mole"))
    {
        FatalError << "Unknown fractionBasis type"
                   << token::SPACE
                   << fractionBasis
                   << nl
                   << "Valid types are: mass or mole."
                   << abort(FatalError);
    }
    fileName chemkinFile(thermo.lookup("CHEMKINFile"));
    fileName thermoFile(thermo.lookup("CHEMKINThermoFile"));
    chemkinReader cr(chemkinFile.expand(), species, thermoFile.expand(), false);
    const HashPtrTable<gasHThermoPhysics>& speciesThermo = cr.speciesThermo(); 
    const HashTable<List<chemkinReader::specieElement> >& specieComposition =
        cr.specieComposition();
    scalarField Cfraction(species.size());
    scalarField Hfraction(species.size());
    scalarField Ofraction(species.size());
    scalarField yFuel(species.size(), 0.0);
    scalarField yOx(species.size(), 0.0);
    scalarField xFuel(species.size(), 0.0);
    scalarField xOx(species.size(), 0.0);
    scalar yC1, yC2;
    scalar yH1, yH2;
    scalar yO1, yO2;
    yC1 = 0;
    yC2 = 0;
    yH1 = 0;
    yH2 = 0;
    yO1 = 0;
    yO2 = 0;
    HashTable<List<chemkinReader::specieElement> >::iterator speciesIter;
    HashPtrTable<gasHThermoPhysics>::iterator speciesThermoIter;
    forAll(species, specieI)
    {
        scalar C = 0;
        scalar H = 0;
        scalar O = 0;
        

        
        forAll(specieComposition[species[specieI]], elemI)
        {
            word name = specieComposition[species[specieI]][elemI].elementName;
            label n = specieComposition[species[specieI]][elemI].nAtoms;
            
            if (name == "C")
            {
                C = atomicWeights["C"]*n/speciesThermo[species[specieI]]->W();
            }
            if (name == "H")
            {
                H = atomicWeights["H"]*n/speciesThermo[species[specieI]]->W();
            }
            if (name == "O")
            {
                O = atomicWeights["O"]*n/speciesThermo[species[specieI]]->W();
            }
        }
        Cfraction[specieI] = C;
        Hfraction[specieI] = H;
        Ofraction[specieI] = O;
        
        if (fractionBasis == "mass")
        {
            if (fuel.found(species[specieI]))
            {
                yFuel[specieI] = readScalar(fuel.lookup(species[specieI]));
            }
            if (oxidizer.found(species[specieI]))
            {
                yOx[specieI] = readScalar(oxidizer.lookup(species[specieI]));
            }
        }
        else
        {
            if (fuel.found(species[specieI]))
            {
                xFuel[specieI] = readScalar(fuel.lookup(species[specieI]));
            }
            if (oxidizer.found(species[specieI]))
            {
                xOx[specieI] = readScalar(oxidizer.lookup(species[specieI]));
            }            
        }
    }
    
    scalar mwf = 0.0;
    scalar mwo = 0.0;
    if (fractionBasis == "mole")
    {
        const scalar mTotf = sum(xFuel);
        const scalar mToto = sum(xOx);
        
        forAll(species, i)
        {
            xFuel[i] /= mTotf;
            xOx[i] /= mToto;
            
            mwf += speciesThermo[species[i]]->W()*xFuel[i];
            mwo += speciesThermo[species[i]]->W()*xOx[i];
        }
        
        forAll(species, i)
        {
            yFuel[i] = xFuel[i]*speciesThermo[species[i]]->W()/mwf;
            yOx[i] = xOx[i]*speciesThermo[species[i]]->W()/mwo;
        }
    }
    forAll(species, specieI)
    {    
        
        
        yC1 += yFuel[specieI]*Cfraction[specieI];
        yH1 += yFuel[specieI]*Hfraction[specieI];
        yO1 += yFuel[specieI]*Ofraction[specieI];

        yC2 += yOx[specieI]*Cfraction[specieI];
        yH2 += yOx[specieI]*Hfraction[specieI];
        yO2 += yOx[specieI]*Ofraction[specieI];
        
    }
     
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Calculating mixture fraction for time "
            << runTime.timeName() 
            << endl;
#       include "createFields.H"
        YC *= 0;
        YH *= 0;
        YO *= 0;
        
         
        forAll(species, specieI)
        {
            YC += Y[specieI]*Cfraction[specieI];
            YH += Y[specieI]*Hfraction[specieI];
            YO += Y[specieI]*Ofraction[specieI];
        }
        Z =
        (
            (
                2.0*(YC - yC2)/atomicWeights["C"]
              + (YH - yH2)/(2.0*atomicWeights["H"])
              - (YO - yO2)/atomicWeights["O"]
            )
           /(
                2.0*(yC1 - yC2)/atomicWeights["C"]
              + (yH1 - yH2)/(2.0*atomicWeights["H"])
              - (yO1 - yO2)/atomicWeights["O"]
            )
        );
        YC.write();
        YH.write();
        YO.write();
        Z.write();
    }
    
    Info<< "\nEnd\n" << endl;
    
    return 0;
}


// ************************************************************************* //
