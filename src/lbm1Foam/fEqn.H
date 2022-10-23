{// BGK Lattice Boltzman equation discretization

	for(unsigned int dI=1; dI<dirNum; dI++)
	{
		for(unsigned int ic=0; ic<internal; ++ic)
		{
			fvScalarMatrix fiEqn
			(
					fvm::ddt(f[dI])
			  + fvm::div(xi[dI], f[dI])
				+ fvm::Sp(omega, f[dI])
				==
					fvc::Sp(omega, feq[dI])
			);

			fiEqn.relax();
			fvConstraints.constrain(fiEqn);
			fiEqn.solve();
			fvConstraints.constrain(f[dI]);
			f[dI].correctBoundaryConditions();

		} // end internal cycles loop

	} // end loop over velocity set

	// solve finally for the resting distribution
	fvScalarMatrix f0Eqn
	(
			fvm::ddt(f[0])
		+ fvm::Sp(omega, f[0])
		==
			fvc::Sp(omega, feq[0])
	);

	f0Eqn.relax();
	fvConstraints.constrain(f0Eqn);
	f0Eqn.solve();
	fvConstraints.constrain(f[0]);
	f[0].correctBoundaryConditions();

	// load particle distributions moments
	rho *= 0.;
	momentum *= 0.;
	forAll(f, dI)
	{
		rho += f[dI];
		momentum += c[dI]*f[dI];
	}

} // end of anonymous namespace
