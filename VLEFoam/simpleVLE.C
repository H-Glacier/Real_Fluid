/*---------------------------------------------------------------------------*\
  =======                |
  \      /  F ield       | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration   |
    \  /    A nd         | Website:  https://openfoam.org
     \/     M anipulation|
\*---------------------------------------------------------------------------*/

#include "simpleVLE.H"
#include "mathematicalConstants.H"
#include "thermodynamicConstants.H"

namespace Foam
{

namespace
{
    inline scalar cubeRoot(const scalar x)
    {
        return x < 0 ? -pow(-x, scalar(1.0/3.0)) : pow(x, scalar(1.0/3.0));
    }

}

simpleVLE::simpleVLE
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
    scalar tpdTol,
    scalar pSatTol,
    label pSatMaxIter
)
:
    thermo_(thermo),
    rho_(rho),
    TcMix_(TcMix),
    PcMix_(PcMix),
    tpd_(tpd),
    pSat_(pSat),
    rhoL_(rhoL),
    rhoV_(rhoV),
    omega_(omega),
    tpdTol_(tpdTol),
    pSatTol_(pSatTol),
    pSatMaxIter_(pSatMaxIter)
{}


scalar simpleVLE::kappa() const
{
    return 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
}


scalar simpleVLE::alpha(const scalar T, const scalar Tc) const
{
    const scalar Tr = max(T/max(Tc, VSMALL), VSMALL);
    return sqr(1.0 + kappa()*(1.0 - sqrt(Tr)));
}


scalar simpleVLE::aCoeff(const scalar Tc, const scalar Pc) const
{
    using constant::thermodynamic::RR;
    return 0.45724*sqr(RR)*sqr(Tc)/max(Pc, VSMALL);
}


scalar simpleVLE::bCoeff(const scalar Tc, const scalar Pc) const
{
    using constant::thermodynamic::RR;
    return 0.07780*RR*Tc/max(Pc, VSMALL);
}


simpleVLE::PRRoots simpleVLE::roots
(
    const scalar p,
    const scalar T,
    const scalar Tc,
    const scalar Pc
) const
{
    using constant::thermodynamic::RR;

    PRRoots r;
    r.n = 0;
    r.z[0] = r.z[1] = r.z[2] = 0;

    if (p <= 0 || T <= 0 || Tc <= 0 || Pc <= 0)
    {
        return r;
    }

    const scalar A = aCoeff(Tc, Pc)*alpha(T, Tc)*p/sqr(RR*T);
    const scalar B = bCoeff(Tc, Pc)*p/(RR*T);

    const scalar a2 = -(1.0 - B);
    const scalar a1 = A - 2.0*B - 3.0*sqr(B);
    const scalar a0 = -(A*B - sqr(B) - pow3(B));

    const scalar Q = (3.0*a1 - sqr(a2))/9.0;
    const scalar Rl = (9.0*a2*a1 - 27.0*a0 - 2.0*pow3(a2))/54.0;
    const scalar D = pow3(Q) + sqr(Rl);

    if (D <= 0 && Q < 0)
    {
        const scalar arg = max(min(Rl/sqrt(-pow3(Q)), 1.0), -1.0);
        const scalar th = acos(arg);
        const scalar qm = 2.0*sqrt(-Q);

        const scalar rootsRaw[3] =
        {
            qm*cos(th/3.0) - a2/3.0,
            qm*cos((th + 2.0*constant::mathematical::pi)/3.0) - a2/3.0,
            qm*cos((th + 4.0*constant::mathematical::pi)/3.0) - a2/3.0
        };

        for (label i=0; i<3; ++i)
        {
            if (rootsRaw[i] > B + SMALL)
            {
                bool duplicate = false;
                for (label j=0; j<r.n; ++j)
                {
                    duplicate = duplicate || mag(rootsRaw[i] - r.z[j]) < 1e-9;
                }
                if (!duplicate)
                {
                    r.z[r.n++] = rootsRaw[i];
                }
            }
        }
    }
    else
    {
        const scalar D05 = sqrt(max(D, scalar(0)));
        const scalar root = cubeRoot(Rl + D05) + cubeRoot(Rl - D05) - a2/3.0;
        if (root > B + SMALL)
        {
            r.z[r.n++] = root;
        }
    }

    for (label i=0; i<r.n; ++i)
    {
        for (label j=i+1; j<r.n; ++j)
        {
            if (r.z[j] < r.z[i])
            {
                const scalar zTmp = r.z[i];
                r.z[i] = r.z[j];
                r.z[j] = zTmp;
            }
        }
    }

    return r;
}


scalar simpleVLE::lnPhi
(
    const scalar Z,
    const scalar p,
    const scalar T,
    const scalar Tc,
    const scalar Pc
) const
{
    using constant::thermodynamic::RR;

    const scalar A = aCoeff(Tc, Pc)*alpha(T, Tc)*p/sqr(RR*T);
    const scalar B = max(bCoeff(Tc, Pc)*p/(RR*T), VSMALL);
    const scalar sqrt2 = sqrt(2.0);

    const scalar zMinusB = max(Z - B, VSMALL);
    const scalar plus = max(Z + (1.0 + sqrt2)*B, VSMALL);
    const scalar minus = max(Z + (1.0 - sqrt2)*B, VSMALL);

    return
        Z - 1.0
      - log(zMinusB)
      - A/(2.0*sqrt2*B)*log(plus/minus);
}


scalar simpleVLE::rhoFromZ
(
    const scalar Z,
    const scalar p,
    const scalar T,
    const scalar W
) const
{
    using constant::thermodynamic::RR;
    return p*W/(max(Z, VSMALL)*RR*max(T, VSMALL));
}


bool simpleVLE::saturationState
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
) const
{
    if (T >= Tc || T <= 0 || Pc <= 0)
    {
        return false;
    }

    scalar pLeft = -1;
    scalar pRight = -1;
    scalar fLeft = 0;
    scalar lastP = -1;
    scalar lastF = 0;
    bool haveLast = false;

    const scalar pMin = max(Pc*1e-8, scalar(1.0));
    const scalar pMax = Pc*(1.0 - 1e-7);

    for (label i=0; i<=160; ++i)
    {
        const scalar s = scalar(i)/160.0;
        const scalar p = pMin*pow(pMax/pMin, s);
        const PRRoots r = roots(p, T, Tc, Pc);

        if (r.n < 2)
        {
            continue;
        }

        const scalar f = lnPhi(r.z[0], p, T, Tc, Pc)
                       - lnPhi(r.z[r.n - 1], p, T, Tc, Pc);

        if (mag(f) < pSatTol_)
        {
            pLeft = pRight = p;
            fLeft = f;
            break;
        }

        if (haveLast && lastF*f <= 0)
        {
            pLeft = lastP;
            pRight = p;
            fLeft = lastF;
            break;
        }

        lastP = p;
        lastF = f;
        haveLast = true;
    }

    if (pLeft < 0)
    {
        return false;
    }

    if (pLeft == pRight)
    {
        ps = pLeft;
    }
    else
    {
        scalar lo = pLeft;
        scalar hi = pRight;
        scalar flo = fLeft;

        for (label iter=0; iter<pSatMaxIter_; ++iter)
        {
            const scalar mid = sqrt(lo*hi);
            const PRRoots r = roots(mid, T, Tc, Pc);
            if (r.n < 2)
            {
                break;
            }

            const scalar fm = lnPhi(r.z[0], mid, T, Tc, Pc)
                            - lnPhi(r.z[r.n - 1], mid, T, Tc, Pc);

            if (mag(fm) < pSatTol_ || mag(log(hi/lo)) < pSatTol_)
            {
                lo = hi = mid;
                break;
            }

            if (flo*fm <= 0)
            {
                hi = mid;
            }
            else
            {
                lo = mid;
                flo = fm;
            }
        }

        ps = sqrt(lo*hi);
    }

    const PRRoots sr = roots(ps, T, Tc, Pc);
    if (sr.n < 2)
    {
        return false;
    }

    zL = sr.z[0];
    zV = sr.z[sr.n - 1];
    rhoL = rhoFromZ(zL, ps, T, W);
    rhoV = rhoFromZ(zV, ps, T, W);

    return rhoL > rhoV && rhoV > 0;
}


label simpleVLE::correct(const volScalarField& e, volScalarField& alphaVap) const
{
    (void)e;

    const scalarField& Tfield = thermo_.T();
    const scalarField& pfield = thermo_.p();
    tmp<volScalarField> tW = thermo_.W();
    const scalarField& Wfield = tW();

    label nTwoPhase = 0;

    forAll(alphaVap, i)
    {
        const scalar T = Tfield[i];
        const scalar p = pfield[i];
        const scalar rhoMix = max(rho_[i], SMALL);
        const scalar Tc = TcMix_[i];
        const scalar Pc = PcMix_[i];
        const scalar W = Wfield[i];

        tpd_[i] = 0.0;
        pSat_[i] = 0.0;
        rhoL_[i] = 0.0;
        rhoV_[i] = 0.0;

        const PRRoots localRoots = roots(p, T, Tc, Pc);
        if (localRoots.n >= 2)
        {
            using constant::thermodynamic::RR;
            const scalar zFeed = p*W/(rhoMix*RR*max(T, VSMALL));
            label feedI = 0;
            scalar best = GREAT;

            for (label ri=0; ri<localRoots.n; ++ri)
            {
                const scalar d = mag(localRoots.z[ri] - zFeed);
                if (d < best)
                {
                    best = d;
                    feedI = ri;
                }
            }

            const scalar lnPhiFeed =
                lnPhi(localRoots.z[feedI], p, T, Tc, Pc);

            scalar minTpd = GREAT;
            for (label ri=0; ri<localRoots.n; ++ri)
            {
                if (ri != feedI)
                {
                    minTpd = min
                    (
                        minTpd,
                        lnPhi(localRoots.z[ri], p, T, Tc, Pc) - lnPhiFeed
                    );
                }
            }
            tpd_[i] = minTpd == GREAT ? 0.0 : minTpd;
        }

        scalar ps = 0;
        scalar rL = 0;
        scalar rV = 0;
        scalar zL = 0;
        scalar zV = 0;

        if (!saturationState(T, Tc, Pc, W, ps, rL, rV, zL, zV))
        {
            alphaVap[i] = (localRoots.n && p < Pc && T < Tc) ? 1.0 : 0.0;
            continue;
        }

        pSat_[i] = ps;
        rhoL_[i] = rL;
        rhoV_[i] = rV;

        if (rhoMix >= rL)
        {
            alphaVap[i] = 0.0;
            continue;
        }

        if (rhoMix <= rV)
        {
            alphaVap[i] = 1.0;
            continue;
        }

        const scalar vMix = 1.0/rhoMix;
        const scalar vL = 1.0/rL;
        const scalar vV = 1.0/rV;

        scalar betaV = (vMix - vL)/max(vV - vL, VSMALL);
        betaV = max(min(betaV, 1.0), 0.0);

        alphaVap[i] = betaV;
        rho_[i] = 1.0/((1.0 - betaV)/rL + betaV/rV);

        if (tpd_[i] < -tpdTol_ || (betaV > SMALL && betaV < 1.0 - SMALL))
        {
            ++nTwoPhase;
        }
    }

    return nTwoPhase;
}

} // End namespace Foam
