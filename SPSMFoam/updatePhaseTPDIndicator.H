// Update phase-separation tendency indicator:
// - root-based component from real-fluid cubic EOS multiple-Z-root detection
// - critical-proximity component to keep sensitivity near dome boundaries
if (mag(tpdRootWeight) > SMALL)
{
    tmp<volScalarField> tForcedRoots = thermo.forcedMaxRootCells();
    forcedMaxRootCells = tForcedRoots();
}
else
{
    forcedMaxRootCells = dimensionedScalar("zero", dimless, 0.0);
}

scalar phaseIndicatorMaxLocal = -GREAT;
forAll(phaseSeparationIndicator, celli)
{
    const scalar forcedRoot = min(max(forcedMaxRootCells[celli], 0.0), 1.0);

    scalar criticalProximity = 0.0;
    if (mag(tpdCriticalWeight) > SMALL)
    {
        const scalar TcSafe = max(TcMix[celli], VSMALL);
        const scalar PcSafe = max(PcMix[celli], VSMALL);
        const scalar thetaT = max((TcSafe - T[celli])/TcSafe, 0.0);
        const scalar thetaP =
            1.0 - min(mag(p[celli] - PcSafe)/(0.05*PcSafe + VSMALL), 1.0);
        criticalProximity = thetaT*thetaP;
    }

    const scalar indicator =
        min
        (
            max
            (
                tpdRootWeight*forcedRoot
              + tpdCriticalWeight*criticalProximity,
                0.0
            ),
            1.0
        );

    phaseSeparationIndicator[celli] = indicator;
    phaseIndicatorMaxLocal = max(phaseIndicatorMaxLocal, indicator);
}

const scalar phaseIndicatorMaxTPD =
    returnReduce(phaseIndicatorMaxLocal, maxOp<scalar>());
(void)phaseIndicatorMaxTPD;
