/*---------------------------------------------------------------------------*\
  =======                |
  \      /  F ield       | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration   |
    \  /    A nd         | Website:  https://openfoam.org
     \/     M anipulation|
\*---------------------------------------------------------------------------*/

// Implementation of the simpleVLE class.  See simpleVLE.H for details.

#include "simpleVLE.H"
// No additional OpenFOAM headers are required beyond those in the header

namespace Foam
{

// Private helper: determine if state is in the two‑phase region.  A cell
// is considered to be within the two‑phase dome if its temperature is
// below the mixture critical temperature and its pressure lies within
// ±5% of the mixture critical pressure.  This simple criterion avoids
// reliance on saturation functions unavailable in this context.
bool simpleVLE::inTwoPhaseRegion(scalar T, scalar p, scalar TcMix, scalar PcMix) const
{
    if (T >= TcMix)
    {
        return false;
    }

    if (mag(p - PcMix) < 0.05 * PcMix)
    {
        return true;
    }

    return false;
}

// Constructor
simpleVLE::simpleVLE(
    const fluidThermo& thermo,
    volScalarField& rho,
    const volScalarField& TcMix,
    const volScalarField& PcMix,
    scalar aL,
    scalar aV,
    scalar Zc
)
:
    thermo_(thermo),
    rho_(rho),
    TcMix_(TcMix),
    PcMix_(PcMix),
    aL_(aL),
    aV_(aV),
    Zc_(Zc)
{}

// Correct the mixture to enforce equilibrium and compute alpha
label simpleVLE::correct(const volScalarField& e, volScalarField& alpha) const
{
    // Retrieve temperature and pressure fields from the thermodynamic
    // package.  The specific internal energy field 'e' is not used in
    // the simplified implementation but retained for interface
    // compatibility.
    const scalarField& Tfield = thermo_.T();
    const scalarField& pfield = thermo_.p();

    label nTwoPhase = 0;

    forAll(alpha, i)
    {
        const scalar Ti    = Tfield[i];
        const scalar pi    = pfield[i];
        const scalar rhoMix = rho_[i];
        const scalar TcMi  = TcMix_[i];
        const scalar PcMi  = PcMix_[i];

        // Determine whether the state is within the two‑phase dome using
        // mixture critical values.  If not, assign alpha based on whether
        // the pressure is below or above the critical pressure and
        // preserve the existing density.
        if (!inTwoPhaseRegion(Ti, pi, TcMi, PcMi))
        {
            if (pi < PcMi)
            {
                // Consider this to be vapour (superheated): alpha=1
                alpha[i] = 1.0;
            }
            else
            {
                // Consider this to be liquid (subcooled): alpha=0
                alpha[i] = 0.0;
            }
            continue;
        }

        // Inside two‑phase region.  Compute approximate saturated
        // densities using mixture critical values.  The critical
        // density is obtained via the definition: rho_c = Pc / (Zc*R*Tc).
        // The gas constant R per unit mass is provided by the
        // thermodynamic package.  A linear law of rectilinear diameters
        // is adopted: rhoL = rho_c*(1 + aL*(1 - T/Tc)),
        // rhoV = rho_c*max(1 - aV*(1 - T/Tc), SMALL), ensuring the
        // densities converge to rho_c at the critical temperature and
        // diverge as temperature decreases.  This simple model does not
        // capture all real‑fluid behaviour but avoids undefined saturation
        // functions and provides a physically plausible lever rule.

        nTwoPhase++;

        const scalar Rgas  = thermo_.R();
        // Guard against division by zero or non‑positive critical values
        scalar TcSafe = max(TcMi, VSMALL);
        scalar PcSafe = max(PcMi, VSMALL);
        scalar rho_c  = PcSafe / (Zc_ * Rgas * TcSafe);
        scalar tau    = 1.0 - Ti / TcSafe;

        // Saturated liquid density
        scalar rhoL = rho_c * (1.0 + aL_ * tau);
        // Saturated vapour density; ensure positivity
        scalar rhoV = rho_c * (1.0 - aV_ * tau);
        if (rhoV < SMALL)
        {
            rhoV = SMALL;
        }

        // Specific volume of the mixture
        const scalar vMix = 1.0 / max(rhoMix, SMALL);

        // Vapour fraction from the lever rule on specific volume
        scalar alphaCell = (vMix - 1.0/rhoL) / ((1.0/rhoV) - (1.0/rhoL));
        // Clamp to [0,1]
        alphaCell = max(min(alphaCell, 1.0), 0.0);

        // Update mixture density using the lever rule.  This ensures
        // consistency between the computed alpha and rho.  Because the
        // saturated densities are computed from the critical values rather
        // than iterative flash, a single update is used rather than
        // iterating to convergence.
        scalar rhoNew = 1.0 / ((1.0 - alphaCell)/rhoL + alphaCell/rhoV);

        alpha[i] = alphaCell;
        rho_[i]  = rhoNew;
    }

    return nTwoPhase;
}

} // End namespace Foam
