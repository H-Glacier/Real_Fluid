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

    // Reference to thermodynamic package
    const fluidThermo& thermo_;

    // Reference to mixture density (input/output)
    volScalarField& rho_;

    // Tolerance and maximum iterations for the UV flash
    scalar tol_;
    label maxIter_;

    // Determine whether a cell is inside the two‑phase dome
    bool inTwoPhaseRegion(scalar T, scalar p) const;

public:

    // Constructors
    simpleVLE(const fluidThermo& thermo, volScalarField& rho,
              scalar tol = 1e-6, label maxIter = 20);

    // Member Functions

    // Perform equilibrium correction.  Given the specific internal energy
    // field e (per unit mass) this routine computes the vapour fraction
    // alpha and updates the mixture density rho_ via the lever rule.  A
    // reference to an alpha field must be supplied; it will be filled
    // with values in [0,1].  When a cell is outside the two‑phase dome
    // (i.e. subcooled liquid or superheated gas) alpha is set to 0 or
    // 1 accordingly and the density is left untouched.  The routine
    // returns the number of cells that underwent two‑phase correction.
    label correct(const volScalarField& e, volScalarField& alpha) const;
};

} // End namespace Foam

#endif
