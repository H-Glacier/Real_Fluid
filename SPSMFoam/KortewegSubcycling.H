// Fixed sub-time-step phase sampling controlled by the TPD indicator.
// A cell is treated as a pure sampled phase during each sub-step.  If the
// local TPD indicator is unstable, the sampled phase flips to the other pure
// phase.  The time-step vapour fraction is the sampled vapour frequency.
{
    initialPhase = alphaVap;

    forAll(tpdSamplePhase, celli)
    {
        tpdSamplePhase[celli] = alphaVap[celli] >= 0.5 ? 1.0 : 0.0;
    }
    tpdSamplePhase.correctBoundaryConditions();

    scalarField vapourSamples(mesh.nCells(), 0.0);

    const dimensionedScalar kappaKortewegCoeff
    (
        "kappaKortewegCoeff",
        dimensionSet(0, 4, -2, 0, 0, 0, 0),
        kappaKorteweg
    );

    const dimensionedScalar deltaTSub
    (
        "deltaTSub",
        runTime.deltaT().dimensions(),
        runTime.deltaTValue()/scalar(nTPDSubCycles)
    );

    UPrime = U;

    label unstableSamplesLocal = 0;
    scalar phaseIndicatorMax = 0.0;

    Info<< "TPD phase sampling subcycles = " << nTPDSubCycles << endl;

    for (label subCycle=0; subCycle<nTPDSubCycles; ++subCycle)
    {
        if (mag(kappaKorteweg) > SMALL)
        {
            gradRho = fvc::grad(rho);
            lapRho = fvc::laplacian(rho);
            KortewegForce = kappaKortewegCoeff*fvc::grad(lapRho);

            UPrime += deltaTSub*KortewegForce
                /max(rho, dimensionedScalar("rhoFloor", rho.dimensions(), SMALL));
        }

        thermo.correct();
        rho = thermo.rho();
        rho.correctBoundaryConditions();

        forAll(TcMix, celli)
        {
            TcMix[celli] = composition.TcMix(celli);
            PcMix[celli] = composition.PcMix(celli);
        }
        forAll(TcMix.boundaryField(), patchi)
        {
            fvPatchScalarField& TcMixPatch = TcMix.boundaryFieldRef()[patchi];
            fvPatchScalarField& PcMixPatch = PcMix.boundaryFieldRef()[patchi];
            forAll(TcMixPatch, facei)
            {
                TcMixPatch[facei] = composition.TcMix(patchi, facei);
                PcMixPatch[facei] = composition.PcMix(patchi, facei);
            }
        }

        {
            #include "updatePhaseTPDIndicator.H"
            phaseIndicatorMax = max(phaseIndicatorMax, phaseIndicatorMaxTPD);
        }

        forAll(tpdSamplePhase, celli)
        {
            if (phaseSeparationIndicator[celli] > phaseTpdThreshold)
            {
                tpdSamplePhase[celli] = 1.0 - tpdSamplePhase[celli];
                ++unstableSamplesLocal;
            }

            vapourSamples[celli] += tpdSamplePhase[celli];
        }
    }

    const scalar invNSub = 1.0/scalar(nTPDSubCycles);
    forAll(alphaVap, celli)
    {
        alphaVap[celli] = min(max(vapourSamples[celli]*invNSub, 0.0), 1.0);
    }

    alphaVap.correctBoundaryConditions();
    tpdSamplePhase.correctBoundaryConditions();

    if (mag(kappaKorteweg) > SMALL)
    {
        U = UPrime;
        U.correctBoundaryConditions();
        splitFallbackActive = 0.0;
    }
    else
    {
        splitFallbackActive = 0.0;
    }

    const label unstableSamples =
        returnReduce(unstableSamplesLocal, sumOp<label>());
    const scalar phaseIndicatorMaxGlobal =
        returnReduce(phaseIndicatorMax, maxOp<scalar>());

    Info<< "TPD phase sampling: max indicator = " << phaseIndicatorMaxGlobal
        << ", unstable cell-substeps = " << unstableSamples
        << ", mean alphaVap = " << alphaVap.weightedAverage(mesh.V()).value()
        << endl;
}
