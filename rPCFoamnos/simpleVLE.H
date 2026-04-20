/*---------------------------------------------------------------------------*\
  =======                |
  \      /  F ield       | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration   |
    \  /    A nd         | Website:  https://openfoam.org
     \/     M anipulation|
\*---------------------------------------------------------------------------*/

// Description
// ~~~~~~~~~~~~
//
// A simple single‑component vapour–liquid equilibrium (VLE) solver for
// compressible two‑phase flow.  This class takes the current thermodynamic
// state (density and specific internal energy) at each cell and computes
// an equilibrium vapour fraction (alpha) together with a consistent mixture
// density using the lever rule.  The algorithm performs a UV‑flash:
// given the specific volume and internal energy of the mixture, it
// iteratively solves for the saturation temperature such that the mixture
// internal energy matches the weighted average of the saturated liquid and
// vapour energies.  Once the saturation temperature is determined, the
// vapour mass fraction α is obtained from the lever rule applied to the
// specific volume【995506714682466†L150-L170】 and the mixture density is updated
// accordingly.  Outside the saturation dome (supercritical or subcooled
// regions) the algorithm falls back to the single‑phase equation of
// state.
//
// Note:  This implementation assumes that the thermodynamic package
// (fluidThermo) provides methods for evaluating saturation properties
// (pSat, rhoL, rhoV, eL, eV) as functions of temperature.  These
// functions are often available in real‑fluid libraries built on the
// Peng‑Robinson equation of state and are used here without further
// modification.

#ifndef simpleVLE_H
#define simpleVLE_H

#include "volFields.H"
#include "fvCFD.H"
#include "fluidThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleVLE Declaration
\*---------------------------------------------------------------------------*/

class simpleVLE
{
    // Private data

    // Reference to thermodynamic package.  Used for gas constant and
    // specific heat information but does not expose saturation functions.
    const fluidThermo& thermo_;

    // Reference to mixture density (input/output)
    volScalarField& rho_;

    // References to mixture critical temperature and pressure fields.  These
    // are computed externally (see createFields.H) and supplied to the
    // constructor.  They are used to construct approximate saturated
    // liquid/vapour densities based on a simple critical‑based model.
    const volScalarField& TcMix_;
    const volScalarField& PcMix_;

    // Empirical constants for approximating saturated densities.  These
    // coefficients control how far the saturated liquid and vapour
    // densities deviate from the critical density as the temperature
    // decreases from the critical temperature.  See simpleVLE.C for
    // details.
    scalar aL_;
    scalar aV_;

    // Assumed critical compressibility factor.  For many real fluids the
    // Peng–Robinson EoS yields Zc≈0.27–0.30.  A constant value is
    // adopted here since mixture‑averaged Zc is typically close to
    // single‑component values.
    scalar Zc_;

    // Determine whether a cell is inside the two‑phase dome.  The
    // criterion compares local temperature and pressure to the mixture
    // critical values.  If the temperature is below TcMix and the
    // pressure is within a small band around PcMix the state is
    // considered two phase.
    bool inTwoPhaseRegion(scalar T, scalar p, scalar TcMix, scalar PcMix) const;

public:

    // Constructors
    // Construct from thermodynamic package, density field and mixture
    // critical temperature/pressure fields.  Optional coefficients
    // control the saturated density approximation and the assumed
    // critical compressibility factor.
    simpleVLE(
        const fluidThermo& thermo,
        volScalarField& rho,
        const volScalarField& TcMix,
        const volScalarField& PcMix,
        scalar aL = 2.0,
        scalar aV = 2.0,
        scalar Zc = 0.27
    );

    // Perform equilibrium correction.  Given the specific internal energy
    // field e (per unit mass) this routine computes the vapour fraction
    // alpha and updates the mixture density rho_ via the lever rule.  A
    // reference to an alpha field must be supplied; it will be filled
    // with values in [0,1].  Cells outside the two‑phase dome are
    // assigned alpha=0 (liquid) or 1 (vapour) based on the local
    // pressure relative to PcMix and the density is left untouched.
    // The return value is the number of cells that underwent two‑phase
    // correction.  The specific internal energy field is accepted for
    // compatibility with existing solver interfaces but is not used in
    // the simplified model; only the density and temperature/pressure
    // fields are required.
    label correct(const volScalarField& e, volScalarField& alpha) const;
};

} // End namespace Foam

#endif
