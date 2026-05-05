
Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

tmp<volScalarField> forced = thermo.forcedMaxRootCells();
label nForced = thermo.nForcedMaxRootCells();

word inertSpecie(thermo.lookup("inertSpecie"));

thermo.correct();

volScalarField& p = thermo.p();
volScalarField& h = thermo.he();
volScalarField& psi = const_cast<volScalarField&>(thermo.psi());
volScalarField& T = thermo.T();

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

volVectorField HbyA ("HbyA", U);


Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K
(
    "K",
    0.5*magSqr(U)
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

// createFields.H
pressureControl pressureLimiter
(
    p,
    rho,
    mesh.solutionDict().subDict("PIMPLE"),
    false
);

const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");

const dimensionedScalar TMin
(
    dimensionedScalar::lookupOrDefault
    (
        "TMin",
        pimpleDict,
        dimTemperature,
        120
    )
);

const dimensionedScalar TMax
(
    dimensionedScalar::lookupOrDefault
    (
        "TMax",
        pimpleDict,
        dimTemperature,
        4500
    )
);

dimensionedScalar v_one("v_one", pow(dimTime,-1.0), 1.0);
volScalarField rAU ("rAU", 1.0/rho/v_one);
rAU.primitiveFieldRef() *= runTime.deltaTValue();
rAU.primitiveFieldRef() /= mesh.V().field();

//- Creates mixture critical temperature field [K]
Info<< "Creating field mixture critical temperature - TcMix \n" << endl;
volScalarField TcMix
(
    IOobject
    (
        "TcMix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TcMix", dimTemperature, 0.0)
);
//- Creates mixture critical pressure field [same units as specie Pc]
Info<< "Creating field mixture critical pressure - PcMix \n" << endl;
volScalarField PcMix
(
    IOobject
    (
        "PcMix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("PcMix", dimPressure, 0.0)
);
forAll(TcMix, celli)
{
    TcMix[celli] = composition.TcMix(celli);
    PcMix[celli] = composition.PcMix(celli);
}
forAll(TcMix.boundaryField(), patchi)
{
    fvPatchScalarField& TcMixp = TcMix.boundaryFieldRef()[patchi];
    fvPatchScalarField& PcMixp = PcMix.boundaryFieldRef()[patchi];
    forAll(TcMixp, facei)
    {
        TcMixp[facei] = composition.TcMix(patchi, facei);
        PcMixp[facei] = composition.PcMix(patchi, facei);
    }
}

volScalarField forcedMaxRootCells
(
    IOobject
    (
        "forcedMaxRootCells",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
volScalarField KChange
(
    "KChange",
    rho*K / runTime.deltaT() * 0.0
);

volScalarField KSource
(
    "KSource",
    p / runTime.deltaT() * 0.0
);

volScalarField c
(
    "cSound",
    sqrt(thermo.Cp()/thermo.Cv() / psi)
);

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

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimTime, 0.0)
);

volVectorField gradp
(
    "gradp",
    fvc::grad(p)
);
volVectorField limGradP 
(
    "limGradP",
    fvc::grad(p)
);

DynamicList<label> badQualityCells(0);

Switch updateEnergyInPISO (false);
word fluxScheme (word::null);
Switch noFluxWeighting(false);

// -------------------------------------------------------------------------
// Create vapour fraction field and instantiate the VLE solver.
// The alpha field represents the vapour mass fraction (between 0 and 1)
// and is initialised to zero.  The simpleVLE solver will update this
// field and the mixture density during the simulation.  The solver
// requires access to the thermodynamic package, the mixture density
// field and the mixture critical temperature/pressure fields, which
// have been created above.  Default coefficient values are used for
// the saturated density approximation.

Info<< "Creating vapour fraction field alphaVap and VLE solver" << endl;

const scalar vleAcentricFactor =
    pimpleDict.lookupOrDefault<scalar>("vleAcentricFactor", 0.0);
const scalar vleTpdTolerance =
    pimpleDict.lookupOrDefault<scalar>("vleTpdTolerance", 1e-8);
const scalar vlePsatTolerance =
    pimpleDict.lookupOrDefault<scalar>("vlePsatTolerance", 1e-8);
const label vlePsatMaxIter =
    pimpleDict.lookupOrDefault<label>("vlePsatMaxIter", 60);

volScalarField alphaVap
(
    IOobject
    (
        "alphaVap",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaVap", dimless, 0.0)
);

volScalarField tpdVLE
(
    IOobject
    (
        "tpdVLE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("tpdVLE", dimless, 0.0)
);

volScalarField pSatVLE
(
    IOobject
    (
        "pSatVLE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("pSatVLE", dimPressure, 0.0)
);

volScalarField rhoLsatVLE
(
    IOobject
    (
        "rhoLsatVLE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoLsatVLE", dimDensity, 0.0)
);

volScalarField rhoVsatVLE
(
    IOobject
    (
        "rhoVsatVLE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoVsatVLE", dimDensity, 0.0)
);

simpleVLE vleSolver
(
    thermo,
    rho,
    TcMix,
    PcMix,
    tpdVLE,
    pSatVLE,
    rhoLsatVLE,
    rhoVsatVLE,
    vleAcentricFactor,
    vleTpdTolerance,
    vlePsatTolerance,
    vlePsatMaxIter
);
