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

#include "laminarTransport.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laminarTransport, 0);
    defineRunTimeSelectionTable(laminarTransport, transportModel);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarTransport::laminarTransport
(
    const volVectorField& U,
    basicMultiComponentMixture& composition,
    moleFraction& moleFraction_,
    psiReactionThermo& thermo,
    const fvMesh& mesh
)
:
    logT
    (
        IOobject
        (
            "logT",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 1.0 )
    ),
    property
    (
        IOobject
        (
            "property",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 1.0 )
    ),
    composition_(composition),
    thermo_(thermo),
    n_(composition.Y().size()),
    Y_(composition.Y()),
    X_(moleFraction_.X()),
    muSpecies_(n_),
    D_(n_*(n_ - 1)/2 + n_),
    mu_
    (
        IOobject
        (
            "mu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, -1, -1, 0, 0, 0, 0), 0 )
    ),
    alpha_
    (
        IOobject
        (
            "alpha",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, -1, -1, 0, 0, 0, 0), 0 )
    ),
    alphaE_
    (
        IOobject
        (
            "alphaE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, -1, -1, 0, 0, 0, 0), 0 )
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0)
    ),    
    V_(n_),
    U_(U),
    a0_(n_),
    a1_(n_),
    a2_(n_),
    a3_(n_),
    d0_(n_),
    d1_(n_),
    d2_(n_),
    d3_(n_),
    transportDict_
    (
        IOdictionary
        (
            IOobject
            (
                "transportProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    ),
    viscousDissipation_(transportDict_.lookupOrDefault("viscousDissipation", false))
{

   forAll(muSpecies_, i)
   {
       const word name = "mu." + Y_[i].name();
       muSpecies_.set
       (
           i,
           new volScalarField
           (
               IOobject
               (
                   name,
                   mesh.time().timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               mesh,
               dimensionedScalar("zero",dimensionSet(1, -1, -1, 0, 0, 0, 0), 0 ) 
               
           )          
           
       );
   }
   
   forAll(Y_, i)
   {
       for (label j = 0; j <= i; j++)
       {
           const label k = i*(i + 1)/2 + j;
           const word Dname = 
               "D." + Y_[i].name() + "." + Y_[j].name();
           D_.set
           (
               k,
               new volScalarField
               (
                  IOobject
                  (
                      Dname,
                      mesh.time().timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::NO_WRITE
                  ),
                  mesh,
                  dimensionedScalar
                  (
                      "zero",
                      dimensionSet(1, -1, -1, 0, 0, 0, 0),
                      0
                  )
               )  
           );
       }
   }   
   forAll(V_, i)
   {
       const word name = "V." + Y_[i].name();
       V_.set
       (
           i,
           new volVectorField
           (
               IOobject
               (
                   name,
                   mesh.time().timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               mesh,
               dimensionedVector
               (
                   "zero",
                   dimensionSet(0, 1, -1, 0, 0, 0, 0),
                   Foam::vector(0, 0, 0)
               )
           )
       );
   }
   read(mesh);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::laminarTransport::read(const fvMesh& mesh)
{
    Info<< "Reading transport properties" << endl;
    
    IOdictionary viscosityPropertiesDict
    (
        IOobject
        (
            
            "viscosityProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    forAll(Y_, i)
    {
        dictionary CoeffsDict = viscosityPropertiesDict.subDict(Y_[i].name());
        a0_[i] = readScalar(CoeffsDict.lookup("a0"));
        a1_[i] = readScalar(CoeffsDict.lookup("a1"));
        a2_[i] = readScalar(CoeffsDict.lookup("a2"));
        a3_[i] = readScalar(CoeffsDict.lookup("a3"));
    }
    
    IOdictionary diffusivityPropertiesDict
    (
        IOobject
        (
            "diffusivityProperties",
             mesh.time().constant(),
             mesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE,
             false           
        )
    );
    forAll(Y_, i)
    {
        dictionary CoeffsDict = diffusivityPropertiesDict.subDict(Y_[i].name());
        forAll(Y_, j)
        {
            dictionary Coeffsij = CoeffsDict.subDict(Y_[j].name());
            d0_[i][j] = readScalar(Coeffsij.lookup("d0"));
            d1_[i][j] = readScalar(Coeffsij.lookup("d1"));
            d2_[i][j] = readScalar(Coeffsij.lookup("d2"));
            d3_[i][j] = readScalar(Coeffsij.lookup("d3"));
        }
    }
}

Foam::tmp<Foam::surfaceScalarField> Foam::laminarTransport::sumJ() const
{
    tmp<surfaceScalarField> tsumJ
    (
        new surfaceScalarField
        (
            IOobject
            (
                "sumJ",
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimensionSet(1, 0, -1, 0, 0, 0, 0),
                0.0
            )
        )    
    );
    
    surfaceScalarField& sumJ = tsumJ();
    
    
    forAll(Y_, specieI)
    {
        sumJ +=
        linearInterpolate
        (
            thermo_.rho()
           *V_[specieI]
           *Y_[specieI]
            
        ) & kappa_.mesh().Sf();
    }
    
    return tsumJ;
}


Foam::tmp<Foam::volSymmTensorField> Foam::laminarTransport::rhoTau() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "rhoTau",
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -mu_*dev(twoSymm(fvc::grad(this->U_)))
        )
    );
}

Foam::tmp<Foam::fvVectorMatrix> Foam::laminarTransport::divRhoTau
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(mu_, U)
      - fvc::div(mu_*dev2(T(fvc::grad(U))))
    );
}



Foam::tmp<volScalarField> Foam::laminarTransport::nu() const
{
    return mu_/thermo_.rho();
}


Foam::tmp<Foam::volScalarField> Foam::laminarTransport::JHs() const
{
    tmp<volVectorField> tJHs
    (
        new volVectorField
        (
            IOobject
            (
                "JHs",
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedVector
            (
                "zero",
                dimensionSet(1, 0, -3, 0, 0, 0, 0),
                vector(0, 0, 0)
            )
        )       
    );
    
    tmp<volScalarField> thSpecie
    (
        new volScalarField
        (
            IOobject
            (
                "hSpecie",
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedScalar("zero", dimEnergy/dimMass, 0.0)
        )
    );
    
    volVectorField& JHs = tJHs();
    volScalarField& hSpecie = thSpecie();
    
    forAll(Y_, specieI)
    {
        forAll(Y_[specieI], cellI)
        {
            const scalar Ti = thermo_.T()[cellI];
            const scalar pi = thermo_.p()[cellI];
            
            hSpecie[cellI] = composition_.Hs(specieI, pi, Ti);
        }
        forAll(thermo_.T().boundaryField(), patchI)
        {
            fvPatchScalarField& pp = thermo_.p().boundaryField()[patchI];
            fvPatchScalarField& pT = thermo_.T().boundaryField()[patchI];
            fvPatchScalarField& ph = hSpecie.boundaryField()[patchI];
            
            forAll(pT, faceI)
            {
                const scalar Ti = pT[faceI];
                const scalar pi = pp[faceI];
                
                ph[faceI] = composition_.Hs(specieI, pi, Ti);
            }
        }
        
        JHs += thermo_.rho()*hSpecie*Y_[specieI]*V_[specieI];
    }
    
    return fvc::div(tJHs);    
}

Foam::tmp<Foam::volScalarField> Foam::laminarTransport::Hconduction() const
{
    tmp<volVectorField> tHconduction
    (
        new volVectorField
        (
            IOobject
            (
                "Hconduction",
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedVector
            (
                "zero",
                dimMass/pow3(dimTime),
                vector(0, 0, 0)
            )
        )
    );
    
    tmp<volScalarField> tH
    (
        new volScalarField
        (
            IOobject
            (
                "H",
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedScalar("zero", dimEnergy/dimMass, 0.0)
        )   
    );
    
    volVectorField& Hconduction = tHconduction();
    volScalarField& H = tH();
    
    forAll(Y_, specieI)
    {
        forAll(Y_[specieI], cellI)
        {
            const scalar Ti = thermo_.T()[cellI];
            const scalar pi = thermo_.p()[cellI];
            
            H[cellI] = composition_.Hs(specieI, pi, Ti);
        }
        forAll(thermo_.T().boundaryField(), patchI)
        {
            fvPatchScalarField& pp = thermo_.p().boundaryField()[patchI];
            fvPatchScalarField& pT = thermo_.T().boundaryField()[patchI];
            fvPatchScalarField& ph = H.boundaryField()[patchI];
            
            forAll(pT, faceI)
            {
                const scalar Ti = pT[faceI];
                const scalar pi = pp[faceI];
                ph[faceI] = composition_.Hs(specieI, pi, Ti);
            }
        }
        
        Hconduction += alpha_*H*fvc::grad(Y_[specieI]);
    }
    
    return fvc::div(Hconduction);
}

Foam::tmp<Foam::volScalarField> Foam::laminarTransport::Econduction() const
{
    tmp<volVectorField> tEconduction
    (
        new volVectorField
        (
            IOobject
            (
                "Econduction",
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedVector
            (
                "zero",
                dimMass/pow3(dimTime),
                vector(0, 0, 0)
            )
        )
    );
    
    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                kappa_.mesh().time().timeName(),
                kappa_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            kappa_.mesh(),
            dimensionedScalar("zero", dimEnergy/dimMass, 0.0)
        )   
    );
    
    volVectorField& Econduction = tEconduction();
    volScalarField& E = tE();
    
    forAll(Y_, specieI)
    {
        forAll(Y_[specieI], cellI)
        {
            const scalar Ti = thermo_.T()[cellI];
            const scalar pi = thermo_.p()[cellI];
            
            E[cellI] = composition_.Es(specieI, pi, Ti);
        }
        forAll(thermo_.T().boundaryField(), patchI)
        {
            fvPatchScalarField& pp = thermo_.p().boundaryField()[patchI];
            fvPatchScalarField& pT = thermo_.T().boundaryField()[patchI];
            fvPatchScalarField& pe = E.boundaryField()[patchI];
            
            forAll(pT, faceI)
            {
                const scalar Ti = pT[faceI];
                const scalar pi = pp[faceI];
                pe[faceI] = composition_.Es(specieI, pi, Ti);
            }
        }
        
        Econduction += alphaE_*E*fvc::grad(Y_[specieI]);
    }
    
    return fvc::div(Econduction);
}

Foam::tmp<Foam::volScalarField> Foam::laminarTransport::W() const
{

    tmp<volScalarField> trW
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("W", Y_[0].group()),
                Y_[0].time().timeName(),
                Y_[0].mesh()
            ),
            Y_[0].mesh(),
            dimensionedScalar("zero", dimless, 0)
        )
    );

    volScalarField& rW = trW();

    forAll(Y_, i)
    {
        rW += Y_[i]/composition_.W(i);
    }

    rW = 1.0/rW;

    return trW;
}

const Foam::label Foam::laminarTransport::index
(
  const label specieI,
  const label specieJ
) const
{
    if(specieI >= specieJ)
    {
        return specieI*(specieI + 1)/2 + specieJ;
    }
    else 
    {
        return specieJ*(specieJ + 1)/2 + specieI;
    }
}

void Foam::laminarTransport::updateBinaryDiffCoeffs()
{
    forAll(Y_, specieI)
    {
        for (label specieJ = 0; specieI >= specieJ; specieJ++)
        {
            const label k = index(specieI, specieJ);
            
            property = 
            (
                d0_[specieI][specieJ]
              + logT
               *(
                    d1_[specieI][specieJ]
                  + logT*(d2_[specieI][specieJ] + d3_[specieI][specieJ]*logT)
                )    
            );
            D_[k] = 
            (
                exp(property)
               *dimensionedScalar
                (
                    "zero",
                    dimensionSet(1, -1, -1, 0, 0, 0, 0),
                    1.0
                )
               /thermo_.p()
               *dimensionedScalar
                (
                    "zero",
                    dimensionSet(1, -1, -2, 0, 0, 0, 0),
                    1.0e+05                    
                )
            );            
        }
    }    
}

Foam::tmp<Foam::volScalarField> Foam::laminarTransport::viscousDissipation() const
{

    tmp<volScalarField> tVD
    (
        new volScalarField
        (
            IOobject
            (
                "viscousDissipation",
                Y_[0].time().timeName(),
                Y_[0].mesh()
            ),
            Y_[0].mesh(),
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0)
        )
    );

    volScalarField& VD = tVD();
    if (viscousDissipation_)
    {
        VD = (rhoTau() && fvc::grad(U_));
    }

    return tVD;
}
// ************************************************************************* //
