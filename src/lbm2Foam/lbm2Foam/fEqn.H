
// solve Lattice Boltzmann Equations

	forAll(f, dI)
	{
		for(unsigned int ic=0; ic<internal; ++ic)
		{
		// dt(f_i) + c_i Grad(f_i) = -1/tau (f_i - f_eq)
			fvScalarMatrix fiEqn
			(
				fvm::ddt(f[dI])
			  + fvm::div(xi[dI], f[dI])
				==
				fvc::Sp(omega, feq[dI])
			  - fvm::Sp(omega, f[dI])
			);

			fiEqn.relax();
			fvConstraints.constrain(fiEqn);
			fiEqn.solve();
			fvConstraints.constrain(f[dI]);

		} // end internal cycles loop

	} // end loop over velocity set

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
