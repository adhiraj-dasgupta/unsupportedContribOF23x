rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        (fvc::interpolate(rho*HbyA) & mesh.Sf())        
      + (
           (pCorr)
          ?(
               rhorAUf*fvc::ddtCorr(rho, U, phi) 
           )
          :(
               0*rhorAUf*fvc::ddtCorr(rho, U, phi)          
           ) 
        )      
    )
);

fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pdEqn
    (
        fvc::ddt(psi) * pReff
      + fvc::div(phiHbyA)
      - fvm::laplacian(rho*rAU, pd)
     ==
        fvOptions(psi, pd, rho.name())
    );

    fvOptions.constrain(pdEqn);

    pdEqn.solve(mesh.solver(pd.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pdEqn.flux();
    }
}
    Info<< "min/max(pd) = "
        << min(pd).value() << ", " << max(pd).value() << endl;

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U = HbyA - rAU*fvc::grad(pd);
U.correctBoundaryConditions();
fvOptions.correct(U);

