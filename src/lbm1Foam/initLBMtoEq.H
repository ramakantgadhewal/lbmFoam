
	Info<< "\nInitialize simulation to equilibrium distribution\n" << endl;

	// start function objects
	runTime.functionObjects().start();

	// Initialize particle mass distribution and
	rho *= 0.;
	forAll(f, dI)
	{
		// initialize as equilibrium distribution
		f[dI] = feq[dI];
		rho += feq[dI];
	}

	// uniform pressure
	p = pRef*dimPres + (rho - density)/ICS2;

	// write IO objects
	forAll(f, dI)
	{
		(f[dI]).write();
	}
	p.write();
