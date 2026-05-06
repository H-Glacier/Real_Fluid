
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
// alphaVap stores the sampled vapour-phase frequency over the TPD subcycles.

Info<< "Creating vapour fraction field alphaVap and TPD phase sampler fields" << endl;

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

// Parameters for Korteweg-stress subcycling
//const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");
const scalar CoNumKortewegMax = pimpleDict.lookupOrDefault<scalar>("CoNumKortewegMax", 0.2);
const label nSubMax = pimpleDict.lookupOrDefault<label>("nKortewegSubMax", 200);
const scalar kappaKorteweg = pimpleDict.lookupOrDefault<scalar>("kappaKorteweg", 0.0);
const scalar phaseTpdThreshold = pimpleDict.lookupOrDefault<scalar>("phaseTpdThreshold", 0.2);
const scalar tpdRootWeight = pimpleDict.lookupOrDefault<scalar>("tpdRootWeight", 1.0);
const scalar tpdCriticalWeight = pimpleDict.lookupOrDefault<scalar>("tpdCriticalWeight", 0.25);
const label nTPDSubCycles = max(pimpleDict.lookupOrDefault<label>("nTPDSubCycles", 10), 1);

volVectorField UPrime
(
    IOobject
    (
        "UPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);

volVectorField KortewegForce
(
    IOobject
    (
        "KortewegForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(p)*0.0
);

volVectorField gradRho
(
    IOobject
    (
        "gradRho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(rho)
);

volScalarField lapRho
(
    IOobject
    (
        "lapRho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::laplacian(rho)
);

volScalarField phaseSeparationIndicator
(
    IOobject
    (
        "phaseSeparationIndicator",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField initialPhase
(
    IOobject
    (
        "initialPhase",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphaVap
);

volScalarField tpdSamplePhase
(
    IOobject
    (
        "tpdSamplePhase",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("tpdSamplePhase", dimless, 0.0)
);

volScalarField splitFallbackActive
(
    IOobject
    (
        "splitFallbackActive",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
