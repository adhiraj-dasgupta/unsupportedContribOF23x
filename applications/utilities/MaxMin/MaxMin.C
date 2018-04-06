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
    MaxMin

Description
    Print the maximum and minimum values of a volumeScalarField.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "pimpleControl.H"
#include "OSspecific.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("fieldName");
    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"
    word fieldName(args.additionalArgs()[0]);
    OFstream outFile(args.path()/fieldName + ".maxmin.out");
    outFile << "Time"
            << token::SPACE
            << "Maximum"
            << token::SPACE
            << "Minimum"
            << endl;
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
           
        if (fieldHeader.headerOk())
        {
          mesh.readUpdate();
          volScalarField field(fieldHeader, mesh);
          outFile << timeDirs[timeI].name()
                  << token::SPACE
                  << max(field).value()
                  << token::SPACE
                  << min(field).value()
                  << endl;          
        }
        else
        {
            Warning << "Field"
                    << token::SPACE
                    << fieldName
                    << token::SPACE
                    << "not found for time"
                    << token::SPACE
                    << timeDirs[timeI].name();
        }

    
     }
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
