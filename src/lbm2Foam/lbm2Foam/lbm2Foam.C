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
    lbm2Foam

Description
    test solver for 2D fluid transport problems employng lattice Boltzman
    equation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  #include "createControl.H"
  #include "createRDeltaT.H"

  #include "createFields.H"
  #include "createFvConstraints.H"
  #include "createFvModels.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // Initilize if user requires
	if(initScheme)
  {
    #include "initLBM.H"
  }
  else
  {
    Info
    << "Simulation initialized to equilibrium distribution.\n"
    << endl;
  }

  Info<< "\nTime loop\n" << endl;

	while(runTime.run())
	{
    runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;

		while (pimple.loop())
		{

			// solve LB equations
			#include "fEqn.H"

		} // end of PIMPLE loop

    // reconstruct macroscopic fileds
    U = momentum/rho;
    p = pRef*dimPres + (rho-density)/ICS2;

    // load equilibrium distribution factor from macroscopic fields
    forAll(feq, dI)
    {
      cDotUBycs2[dI] = ICS2*(c[dI] & U);
      uEqFactor[dI] = W[dI]*( 1.0
                            + cDotUBycs2[dI]*(1. + 0.5*cDotUBycs2[dI])
                            - 0.5*ICS2*(U&U)
                            );
      feq[dI] = rho*uEqFactor[dI];
    }

    runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "    ClockTime = " << runTime.elapsedClockTime() << " s" << nl << nl
        << "System mass = "
        << (fvc::domainIntegrate(rho)).value() << " kg" << nl
        << endl;

	} // end of time loop

	Info<< "End.\n" << endl;


  return 0;
} // end of main


// ************************************************************************* //
