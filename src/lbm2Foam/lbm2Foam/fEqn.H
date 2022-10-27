{// BGK Lattice Boltzman equation discretization

	for(label ic=0; ic<internal; ++ic)
	{
		for(label dI=1; dI<dirNum; dI++)
		{
			fvScalarMatrix fiEqn
			(
					fvm::ddt(f[dI])
				+ fvm::div(xi[dI], f[dI])
				+ omega*fvm::Sp(1.-uEqFactor[dI], f[dI])
				==
					omega*fvc::Sp(1.-uEqFactor[dI], f[dI])
				- fvc::Sp(omega, f[dI])
				+	fvc::Sp(omega, feq[dI])
			);

			fiEqn.relax();
			fvConstraints.constrain(fiEqn);
			fiEqn.solve();
			fvConstraints.constrain(f[dI]);

			//load distribution gradient
			gradf[dI] = fvc::grad(f[dI]);

		} // end loop over velocity set

	} // end internal cycles loop

	// solve finally for the resting distribution
	fvScalarMatrix f0Eqn
	(
				fvm::ddt(f[0])
			+ omega*fvm::Sp(1.-uEqFactor[0], f[0])
			==
				omega*fvc::Sp(1.-uEqFactor[0], f[0])
			- fvc::Sp(omega, f[0])
			+	fvc::Sp(omega, feq[0])
	);

	f0Eqn.relax();
	fvConstraints.constrain(f0Eqn);
	f0Eqn.solve();
	fvConstraints.constrain(f[0]);

	// load distribution gradient
	gradf[0] = fvc::grad(f[0]);

	// load particle distributions moments
	rho *= 0.;
	momentum *= 0.;
	forAll(f, dI)
	{
		rho += f[dI];
		momentum += c[dI]*f[dI];
	}

	/* load equilibrium distributions:
		explicit treatment of the velocity field
		-> uEqFactor taken from the prevoius time step
	*/
	forAll(feq, dI)
	{
		feq[dI] = uEqFactor[dI]*rho;
	}

} // end of anonymous namespace
