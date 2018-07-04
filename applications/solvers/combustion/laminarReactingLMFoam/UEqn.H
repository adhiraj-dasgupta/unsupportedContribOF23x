     fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + transport.divRhoTau(U)
     ==
        rho*g
      + fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(pd));

        fvOptions.correct(U);
    }
