Info<< "Creating reaction model\n" << endl;

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

word inertSpecie(thermo.lookup("inertSpecie"));

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


IOdictionary chemistryProperties
(
    IOobject
    (
        "chemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
dimensionedScalar pReff(chemistryProperties.lookup("pReff"));
volScalarField& p = thermo.p();
p = pReff;
const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();

#include "compressibleCreatePhi.H"

Info<< "Creating field pd\n" << endl;
volScalarField pd
(
    IOobject
    (
        "pd",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the reaction model
reaction->setTurbulence(turbulence());


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
    fields.add(X[i]);
}
fields.add(thermo.he());

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

autoPtr<laminarTransport> ptransport
(
  laminarTransport::New
  (
      U,
      composition,
      moleFraction_,
      thermo,
      mesh
  )
);
laminarTransport& transport = ptransport();

Switch pCorr(thermo.lookupOrDefault("pCorr", false));
if (pCorr)
{
    Info<< "Extra term added to flux in p equation" << endl;
}

PtrList<volScalarField> ROP_(Y.size());
forAll(Y, i)
{
    const word RRname = "ROP." + Y[i].name();
    ROP_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                RRname,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );    
}
