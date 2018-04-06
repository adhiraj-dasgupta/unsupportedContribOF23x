{
    volScalarField& he = thermo.he();
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he, "div(phi,he)")
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? (
              fvc::div
              (
                  fvc::absolute(phi/fvc::interpolate(rho), U),
                  p,
                  "div(phiv,p)"
              )
            - fvm::laplacian(transport.alphaE(), he)
            + transport.Econduction() 
            )
          : (
                -dpdt
               - fvm::laplacian(transport.alpha(), he)
               + transport.Hconduction()
            )
        )
      + transport.JHs()      
     ==
        rho*(U&g)  
      + reaction->Sh()
      + (U & fvc::div(transport.rhoTau()))
      + transport.viscousDissipation()
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
