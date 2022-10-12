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

	  #include "initLBM.H"

    Info<< "\nTime loop\n" << endl;

	while(runTime.run())
	{
        runTime++;

		while (pimple.loop())
		{
			runTime++;

			Info<< "Time = " << runTime.timeName() << nl << endl;

			// solve LB equations
			#include "fEqn.H"

			// load macroscopic fields from particle distributions
			rho *= 0.;
			momentum *= 0.;
			forAll(f, dI)
			{
				rho += f[dI];
				momentum += c[dI]*f[dI];
			}
			U = momentum/rho;
			p = pRef*dimPres + (rho-density)/ICS2;

			// load equilibrium distributions from macroscopic fields
			forAll(feq, dI)
			{
				feq[dI] = W[dI]*rho*( 1.0
									+ ICS2*(c[dI]&U)
									- 0.5*ICS2*(U&U)
									+ 0.5*ICS4*(c[dI]&U)*(c[dI]&U)
									);
			}

			runTime.write();
		}

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "    ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

	}

	Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
