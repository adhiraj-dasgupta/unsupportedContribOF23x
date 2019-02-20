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
    conditionalMean

Description
    Given two volScalarFields, calculate the conditional mean of the second one
    with respect to the first using the given number of bins.
    Usage:
    conditionalMean <field1> <field2> <nbins>

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "bound.H"
#include "pimpleControl.H"
#include "makeGraph.H"
#include "interpolateXY.H"
#include "OSspecific.H"
#include <cstdlib>



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("fieldName0");
    argList::validArgs.append("fieldName1");
    argList::validArgs.append("nBins");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    word fieldName0(args.additionalArgs()[0]);
    word fieldName1(args.additionalArgs()[1]);
    word bins(args.additionalArgs()[2]);
    char* pEnd;
    scalar nBins = strtod(bins.c_str(), &pEnd);
    scalarField var0(nBins);
    scalarField var1(nBins);
    
    const word& gFormat = runTime.graphFormat();
    Switch isVar0 = false;
    Switch isVar1 = false;
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader0
        (
            fieldName0,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
           
        if (fieldHeader0.headerOk())
        {
          mesh.readUpdate();
          isVar0 = true;
        }
        else
        {
            Warning << "Field"
                    << token::SPACE
                    << fieldName0
                    << token::SPACE
                    << "not found for time"
                    << token::SPACE
                    << timeDirs[timeI].name();
          isVar0 = false;
          continue;
        }

        IOobject fieldHeader1
        (
            fieldName1,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
           
        if (fieldHeader1.headerOk())
        {
          mesh.readUpdate();
          isVar1 = true;
          volScalarField field1(fieldHeader1, mesh);
        }
        else
        {
            Warning << "Field"
                    << token::SPACE
                    << fieldName1
                    << token::SPACE
                    << "not found for time"
                    << token::SPACE
                    << timeDirs[timeI].name();
            isVar1 = false;
            continue;
        }
        
        volScalarField field0(fieldHeader0, mesh);
        volScalarField field1(fieldHeader1, mesh);
        labelField nItems(nBins, 0);
        const scalar max0(max(field0).value());
        const scalar min0(min(field0).value());
        const scalar width = (max0 - min0)/nBins;
        forAll(var0, i)
        {
            const scalar vleft = min0 + i*width;
            const scalar vright = min0 + (i+1)*width;
            var0[i] = (vleft + vright)/2.0;
        }
        var1 = 0.0;
        forAll(field0, cellI)
        {
            const scalar f0 = field0[cellI];
            const scalar f1 = field1[cellI];
            
            for (label i = 0; i < nBins; i++)
            {
                if (
                       (f0 > min0 + i*width)
                    && (f0 <= min0 + (i + 1)*width)
                   )
                {
                   nItems[i]++;
                   var1[i] += f1;
                   continue;
                }                  
            }
        }
        forAll(var1, i)
        {
//            var1[i] /= (nItems[i] + SMALL);
              var1[i] = var1[i]/(nItems[i] + SMALL);
        }
        
        fileName path(field0.rootPath()/field0.caseName()/"conditionalMean"/field0.instance());
        mkDir(path);
        
        const word name = fieldName0 + "." + fieldName1;
        makeGraph(var0, var1, name, path, gFormat);
        
        
    
     }
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
