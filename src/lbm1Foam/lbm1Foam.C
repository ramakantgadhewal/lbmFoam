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
    lbm1Foam

Description
    test solver for 1D transport problems

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"

    #include "createFields.H"
    #include "createFvConstraints.H"
    #include "createFvModels.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  if (initScheme == 0)
  {
    #include "initLBMtoEq.H"
  }
  else
	{
    #include "initLBM.H"
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

      Info<< "System mass = "
          << (fvc::domainIntegrate(rho)).value() << nl << endl;

      // reconstruct macroscopic fileds
      U = momentum/rho;
      p = pRef*dimPres + (rho-density)/ICS2;

      const dimensionedVector gravityDU
      (
        dimVelocity,
        (9.81*deltaT)*Foam::vector(-1,0,0)
      );

      // load equilibrium distributions from macroscopic fields
      forAll(feq, dI)
      {
        cDotU[dI] 	= (c[dI] & (U + gravityDU));
        uEqFactor[dI] = W[dI]*( 1.0
                                    + ICS2*cDotU[dI]
                                    - 0.5*ICS2*(U&U)
                                    + 0.5*ICS4*cDotU[dI]*cDotU[dI]
                                  );
        feq[dI] 		= rho*uEqFactor[dI];
      }

		} // end of PIMPLE loop

    runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "    ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

	}

	Info<< "End\n" << endl;

  return 0;
}


// ************************************************************************* //
