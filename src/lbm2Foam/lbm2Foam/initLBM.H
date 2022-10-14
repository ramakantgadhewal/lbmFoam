
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Info<< "\nInitialize thermo-hydraulic problem.\n" << endl;

	// start function objects
	runTime.functionObjects().start();

	// compute domain volume
	scalar volume =
	(
		fvc::domainIntegrate(rho/rho)
	).value();
	Info<< "Volume = " << volume << nl << endl;

	// Initialize particle mass distribution and
	rho*=0.;
	momentum*=0.;
	scalar uIn=
		runTime.controlDict().lookupOrDefault<scalar>("initUfactor", 1);
	forAll(f, dI)
	{
		f[dI] = W[dI]*rho_prev*( 1.0
							  + ICS2*(c[dI]&U)*uIn
							  + 0.5*ICS4*(c[dI]&U)*uIn*(c[dI]&U)*uIn
							  - 0.5*ICS2*(U&U)*uIn
							  );
		rho += f[dI];
		momentum += c[dI]*f[dI];
	}
	scalar errIC = 0.;

	Info<< "Starting density error = "
			<< average(mag(1.-rho/rho_prev)).value() << nl << endl;

	// loop to fit density field and velocity field
	scalar tolIC =
		runTime.controlDict().lookupOrDefault<scalar>("initTol", 1e-9);
	unsigned int iterICmax =
		runTime.controlDict().lookupOrDefault<scalar>("initIterMax", 1e4);
	unsigned int iIC = 0;
	while(errIC>tolIC && iIC<iterICmax)
	{
		// load previous density
		rho_prev = rho;

		// equilibrium ditribution
		forAll(feq, dI)
		{
			feq[dI] = W[dI]*rho*( 1.0
								+ ICS2*(c[dI]&U)
								+ 0.5*ICS4*(c[dI]&U)*(c[dI]&U)
								- 0.5*ICS2*(U&U)
								);
		}

		// solve LBE
		#include "fEqn.H"

		// evaluate new density and velocity
		rho*=0.;
		momentum*=0.;
		forAll(f, dI)
		{
			rho+=f[dI];
			momentum+=c[dI]*f[dI];
		}

		// compute relative incremental error
		errIC =
		(
			average(1.-mag(momentum)/rho/mag(U)) // velocity error
			+ average(1.-rho/rho_prev) // pressure error
		).value();

		iIC++;
	}

	// evaluate initial pressure distribution
	p = pRef*dimPres + (rho-density)/ICS2;


/*	Info<< "Initialization iterations = " << iIC << nl
			<< "Velocity error reached = " << errIC << nl
			<< "Initial pressure error = " << average(mag(1.0-rho/rho_prev)).value() << nl
			<< "Mass loss = " << fvc::domainIntegrate(density-rho).value() << nl
			<< "Mass fraction = " << average(rho/rho0).value() << nl
			<< "Velocity fraction = " << average(mag(momentum)/rho/mag(U)).value() << nl << endl;
*/
	// overwrite discontinuous velocity
	//U = momentum/rho;
	//U.write();

	// write IO
	forAll(f, dI)
	{
		(f[dI]).write();
	}
	p.write();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Info<< "\nInitialize neutronic problem.\n" << endl;

	// equilibrium distributions
	forAll(feq_flux, dI)
	{
		feq_flux[dI] 	= W[dI]*flux;
		f_flux[dI] 		= W[dI]*flux;
		(f_flux[dI]).write();
	}

	// set precursors to steady state
	delSource *= 0.;
	forAll(prec, pI)
	{
		prec[pI] = beta[pI]/lam[pI]*NSF0_0*flux;
		delSource += cl*lam[pI]*prec[pI];
		prec[pI].write();
	}
