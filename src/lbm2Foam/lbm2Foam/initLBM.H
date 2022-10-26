
	// initialization options
	const scalar tolIC =
		runTime.controlDict().lookupOrDefault<scalar>("initTol", 1e-12);

	const unsigned int iterICmax =
		runTime.controlDict().lookupOrDefault<scalar>("initIterMax", 1e3);

	const unsigned int iterICerrEval =
		runTime.controlDict().lookupOrDefault<scalar>("initErrEval", 1e1);

	Info<< "\nInitialize simulation\n" << endl;

	// start function objects
	runTime.functionObjects().start();

	// Initialize particle mass distribution and
	rho*=0.;
	forAll(f, dI)
	{
		// initialize as zero velocity - uniform pressure
		f[dI] = W[dI]*rho_prev;
		rho += f[dI];
	}

	// loop to fit density field and velocity field
	scalar errIC = 1;
	unsigned int iIC = 0;

	while(errIC>tolIC && iIC<iterICmax)
	{
		// evaluate equilibrium distribution at fixed
		// velocity field U
		forAll(feq, dI)
		{
			feq[dI] = uEqFactor[dI]*rho;
		}

		// solve LBE
		#include "fEqn.H"

		// compute relative incremental error
		if(iIC%iterICerrEval==0)
		{
			errIC =
			(
				average(mag(1.-rho/rho_prev)) // pressure error
			).value();

			// load previous density
			rho_prev = rho;
		}
		iIC++;
	}

	Info<< " Initialization iterations = "
			<< iIC << nl
			<< " Initializzation pressure errror = "
			<< errIC << nl << endl;

	// evaluate initial pressure distribution
	p = pRef*dimPres + (rho-density)/ICS2;

	// write IO objects
	forAll(f, dI)
	{
		(f[dI]).write();
	}
	p.write();
