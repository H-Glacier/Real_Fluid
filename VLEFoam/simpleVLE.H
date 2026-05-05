/*---------------------------------------------------------------------------*\
  =======                |
  \      /  F ield       | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration   |
    \  /    A nd         | Website:  https://openfoam.org
     \/     M anipulation|
\*---------------------------------------------------------------------------*/

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
    const fluidThermo& thermo_;
    volScalarField& rho_;

    const volScalarField& TcMix_;
    const volScalarField& PcMix_;

    volScalarField& tpd_;
    volScalarField& pSat_;
    volScalarField& rhoL_;
    volScalarField& rhoV_;

    scalar omega_;
    scalar tpdTol_;
    scalar pSatTol_;
    label pSatMaxIter_;

    struct PRRoots
    {
        label n;
        scalar z[3];
    };

    scalar kappa() const;
    scalar alpha(const scalar T, const scalar Tc) const;
    scalar aCoeff(const scalar Tc, const scalar Pc) const;
    scalar bCoeff(const scalar Tc, const scalar Pc) const;

    PRRoots roots
    (
        const scalar p,
        const scalar T,
        const scalar Tc,
        const scalar Pc
    ) const;

    scalar lnPhi
    (
        const scalar Z,
        const scalar p,
        const scalar T,
        const scalar Tc,
        const scalar Pc
    ) const;

    scalar rhoFromZ
    (
        const scalar Z,
        const scalar p,
        const scalar T,
        const scalar W
    ) const;

    bool saturationState
    (
        const scalar T,
        const scalar Tc,
        const scalar Pc,
        const scalar W,
        scalar& ps,
        scalar& rhoL,
        scalar& rhoV,
        scalar& zL,
        scalar& zV
    ) const;

public:

    simpleVLE
    (
        const fluidThermo& thermo,
        volScalarField& rho,
        const volScalarField& TcMix,
        const volScalarField& PcMix,
        volScalarField& tpd,
        volScalarField& pSat,
        volScalarField& rhoL,
        volScalarField& rhoV,
        scalar omega,
        scalar tpdTol = 1e-8,
        scalar pSatTol = 1e-8,
        label pSatMaxIter = 60
    );

    // Single-component Peng-Robinson VLE correction.  The pure-fluid TPD is
    // min(ln(phiTrial/phiFeed)) over the alternate PR root at the local T,p.
    // Saturated liquid/vapour densities come from phiL(T,pSat)=phiV(T,pSat).
    label correct(const volScalarField& e, volScalarField& alpha) const;
};

} // End namespace Foam

#endif
