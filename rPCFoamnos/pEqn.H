rAU = 1.0 / UEqn.A();
HbyA = UEqn.H() * rAU;
HbyA.boundaryFieldRef() == U.boundaryField();

#include "updateCentralWeights.H"
#include "psiUpdateCentralFields.H"

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn_own
    (
        fvm::div(phid_own,p) - fvm::laplacian(Dp_own,p)
    );
    
    fvScalarMatrix pEqn_nei
    (
        fvm::div(phid_nei,p) - fvm::laplacian(Dp_nei,p)
    );
    
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi,p)
        +
        pEqn_own
        +
        pEqn_nei
        ==
        fvOptions(psi,p,rho.name())
    );
    
    fvOptions.constrain(pEqn);
    
    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
    p.max(dimensionedScalar("pMin", dimPressure, 8e6));
    p.min(dimensionedScalar("pMax", dimPressure, 2.5e7));
  
    
    if (pimple.finalNonOrthogonalIter())
    {
        phi_own = pEqn_own.flux();
        phi_nei = pEqn_nei.flux();
        phi = phi_own + phi_nei;
    }
}

p_own = fvc::interpolate(p, own, "reconstruct(p)");
p_nei = fvc::interpolate(p, nei, "reconstruct(p)");

rho_own = p_own*psi_own;
rho_nei = p_nei*psi_nei;
correctCentralACMIInterpolation(rho_nei);

gradp = fvc::div((alpha_own*p_own + alpha_nei*p_nei)*mesh.Sf());
#include "limitPressureGradient.H"

U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);

Info << "max(U): " << max(U).value() << endl;
Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

#include "massEqn.H"

volScalarField rhoT = thermo.rho();

#include "centralContinuityErrs.H"

rho = rhoT;

        // Perform vapour–liquid equilibrium correction on the updated
        // thermodynamic fields.  This routine computes the vapour
        // fraction alpha and updates the mixture density via the lever
        // rule using the simplified VLE model.  The specific internal
        // energy is retrieved from the thermodynamic package but is not
        // used internally by the simpleVLE implementation; it is passed
        // for interface compatibility.  See simpleVLE.C for details.
        {
            tmp<volScalarField> tE = thermo.e();
            // Pass the vapour fraction field alphaVap to the VLE solver
            vleSolver.correct(tE(), alphaVap);
        }


