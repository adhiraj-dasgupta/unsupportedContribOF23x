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

Class
    Foam::QSS

Description
    A solver for chemistry when the mechanism is available as a Fortran ckwyp
    routine. Reduced mechanisms like these contain QSS species which are not
    transported, but only occur in the chemical source terms. Since fewer
    species are transported, this can lead to a saving in computational
    effort.

SourceFiles
    QSS.C

\*---------------------------------------------------------------------------*/

#ifndef QSS_H
#define QSS_H

#include "chemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

        //- the external routine to compute the reaction rates
        extern "C"
        {
            void ckwyp_
            (
                double* p,
                double* T,
                double  [],
                int* ickwrk,
                double* rckwrk,
                double [] 
            );
        }

/*---------------------------------------------------------------------------*\
                       Class QSS Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class QSS
:
    public chemistrySolver<ChemistryModel>
{
    // Private data

        dictionary coeffsDict_;
        autoPtr<ODESolver> odeSolver_;

        // Solver data
        mutable scalarField cTp_;
        mutable double* Ytmp0;
        mutable double* dYdt;


public:

    //- Runtime type information
    TypeName("QSS");


    // Constructors

        //- Construct from mesh
        QSS(const fvMesh& mesh);


    //- Destructor
    virtual ~QSS();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;



        void derivatives
        (
            const scalar t,
            const scalarField& c,
            scalarField& dcdt
        ) const;

        void jacobian
        (
            const scalar t,
            const scalarField& c,
            scalarField& dcdt,
            scalarSquareMatrix& dfdc
        ) const;
        
        tmp<scalarField> omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;
        
        template<class ThermoType>
        scalar omega
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
        ) const;

        virtual scalar omegaI
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
        ) const;
        
        virtual tmp<DimensionedField<scalar, volMesh> > calculateRR
        (
            const label reactionI,
            const label specieI    
        ) const;
        
        virtual tmp<volScalarField> tc() const;
        
        //- Calculates the reaction rates
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "QSS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
