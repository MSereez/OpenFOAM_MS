/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    transientSimpleFoam

Description
   Large Time Step transient Solver for incompressible flows using SIMPLE algorithm

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include <stdlib.h>     /* strtod */
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"


    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

       Info<< "deltaT  = " <<  runTime.deltaTValue() << endl;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        // --- Pressure-velocity PIMPLE corrector loop
         while (piso.correct())
        {
        	// Solve the Momentum equation U

        	tmp<fvVectorMatrix> tUEqn
        	(
        	    fvm::ddt(U) + fvm::div(phi, U)
        	  + turbulence->divDevReff(U)
        	);

        	fvVectorMatrix& UEqn = tUEqn.ref();
        	UEqn.relax();

        	solve(UEqn == -fvc::grad(p));


         //   p.boundaryFieldRef() == p.boundaryField();

            p.boundaryFieldRef().updateCoeffs();

            volScalarField rAU(1.0/UEqn.A());
            U = rAU * UEqn.H();
            tUEqn.clear();
            phi = fvc::interpolate(U) & mesh.Sf();
            adjustPhi(phi, U, p);

            // Store pressure for under-relaxation
            p.storePrevIter();

            tmp<volScalarField> rAtU(rAU);

                while (piso.correctNonOrthogonal())
                {

            	    // Pressure corrector
            	    fvScalarMatrix pEqn
            	    (
            	        fvm::laplacian(rAtU(), p) == fvc::div(phi)
            	    );

            	    pEqn.setReference(pRefCell, pRefValue);

            	    pEqn.solve();

            	    if (piso.finalNonOrthogonalIter())
            	    {
            	        phi -= pEqn.flux();
            	    }
            	}

            	#include "continuityErrs.H"

            	// Explicitly relax pressure for momentum corrector
            	p.relax();

                // Momentum corrector
                U -= rAtU()*fvc::grad(p);
                U.correctBoundaryConditions();
        }

            laminarTransport.correct();
            turbulence->correct();


            runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
      }

    Info<< "End\n" << endl;

    return 0;
 }



// ************************************************************************* //
