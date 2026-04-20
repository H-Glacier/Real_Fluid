/*---------------------------------------------------------------------------*\
  =======                |
  \      /  F ield       | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration   |
    \  /    A nd         | Website:  https://openfoam.org
     \/     M anipulation|
\*---------------------------------------------------------------------------*/

// Implementation of the simpleVLE class.  See simpleVLE.H for details.

#include "simpleVLE.H"
#include "thermophysicalProperties.H"

namespace Foam
{

// Private helper: determine if state is in two‑phase region
bool simpleVLE::inTwoPhaseRegion(scalar T, scalar p) const
{
    // Use a crude criterion: compare the local pressure to the saturation
    // pressure at the same temperature.  If the pressure is within
    // ±5% of the saturation pressure and below the critical temperature
    // then treat the cell as in the two‑phase region.  This heuristic is
    // commonly used when detailed spinodal detection is not available【995506714682466†L150-L170】.
    scalar pSat = thermo_.pSat(T);
    scalar Tc   = thermo_.Tc();

    if (T >= Tc)
    {
        return false;
    }

    if (mag(p - pSat) < 0.05 * pSat)
    {
        return true;
    }

    return false;
}

// Constructor
simpleVLE::simpleVLE(const fluidThermo& thermo, volScalarField& rho,
                     scalar tol, label maxIter)
:
    thermo_(thermo),
    rho_(rho),
    tol_(tol),
    maxIter_(maxIter)
{}

// Correct the mixture to enforce equilibrium and compute alpha
label simpleVLE::correct(const volScalarField& e, volScalarField& alpha) const
{
    const scalarField& T = thermo_.T();
    const scalarField& p = thermo_.p();

    label nTwoPhase = 0;

    forAll(alpha, i)
    {
        scalar Ti = T[i];
        scalar pi = p[i];
        scalar eMix = e[i];
        scalar rhoMix = rho_[i];

        // Determine whether the state is within the two‑phase dome
        if (!inTwoPhaseRegion(Ti, pi))
        {
            // Single‑phase: set alpha based on whether pressure is above or below
            // saturation line.  If p < pSat, the mixture is vapour; else liquid.
            scalar pSat = thermo_.pSat(Ti);
            if (pi < pSat)
            {
                alpha[i] = 1.0; // vapour
            }
            else
            {
                alpha[i] = 0.0; // liquid
            }
            // Keep the current mixture density
            continue;
        }

        // Otherwise perform UV flash to find equilibrium T and alpha
        nTwoPhase++;

        // Initial guess for saturation temperature is the current temperature
        scalar Tguess = Ti;

        // Specific volume of the mixture
        const scalar vMix = 1.0 / rhoMix;

        for (label iter = 0; iter < maxIter_; ++iter)
        {
            // Evaluate saturated properties at Tguess
            scalar pSat = thermo_.pSat(Tguess);
            scalar rhoL = thermo_.rhoL(Tguess);
            scalar rhoV = thermo_.rhoV(Tguess);
            scalar eL   = thermo_.eL(Tguess);
            scalar eV   = thermo_.eV(Tguess);

            // Compute vapour mass fraction from lever rule on specific volume
            scalar alphaGuess = (vMix - 1.0/rhoL) / ((1.0/rhoV) - (1.0/rhoL));
            alphaGuess = max(min(alphaGuess, 1.0), 0.0);

            // Predicted internal energy
            scalar ePred = (1.0 - alphaGuess)*eL + alphaGuess*eV;

            // Residual: difference between predicted and actual energy
            scalar residual = ePred - eMix;

            if (mag(residual) < tol_)
            {
                // Converged
                break;
            }

            // Finite difference to approximate derivative d(ePred)/dT
            scalar dT = 1.0; // one Kelvin increment
            scalar Tfd = Tguess + dT;

            scalar rhoLfd = thermo_.rhoL(Tfd);
            scalar rhoVfd = thermo_.rhoV(Tfd);
            scalar eLfd   = thermo_.eL(Tfd);
            scalar eVfd   = thermo_.eV(Tfd);
            scalar alphaFd = (vMix - 1.0/rhoLfd) / ((1.0/rhoVfd) - (1.0/rhoLfd));
            alphaFd = max(min(alphaFd, 1.0), 0.0);
            scalar ePredFd = (1.0 - alphaFd)*eLfd + alphaFd*eVfd;
            scalar dE_dT = (ePredFd - ePred) / dT;

            // Prevent division by zero
            if (mag(dE_dT) < SMALL)
            {
                break;
            }

            // Newton update on Tguess
            Tguess -= residual / dE_dT;

            // Clamp the temperature to a reasonable range to avoid divergence
            scalar Tc = thermo_.Tc();
            scalar Tmin = thermo_.TMin();
            scalar Tmax = Tc;
            Tguess = min(max(Tguess, Tmin), Tmax);
        }

        // Evaluate final saturated properties
        scalar rhoL = thermo_.rhoL(Tguess);
        scalar rhoV = thermo_.rhoV(Tguess);
        scalar alphaFinal = (vMix - 1.0/rhoL) / ((1.0/rhoV) - (1.0/rhoL));
        alphaFinal = max(min(alphaFinal, 1.0), 0.0);

        // Update mixture density using the lever rule
        scalar rhoNew = 1.0 / ((1.0 - alphaFinal)/rhoL + alphaFinal/rhoV);

        // Store results
        alpha[i] = alphaFinal;
        rho_[i] = rhoNew;
    }

    return nTwoPhase;
}

} // End namespace Foam
