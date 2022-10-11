// Lattice Boltzman equation discretization

	forAll(f, dI)
	{
		for(int ic=0; ic<internal; ++ic)
		{
			fvScalarMatrix fiEqn
			(
				fvm::ddt(f[dI])
			  + fvm::div(xi[dI], f[dI])
				==
				fvc::Sp(omega, feq[dI])
			  - fvm::Sp(omega, f[dI])
			//+  fvOptions(f[dI])
			);

			fiEqn.relax();
			fvOptions.constrain(fiEqn);
			fiEqn.solve();
			fvOptions.correct(f[dI]);

		} // end internal cycles loop

	} // end loop over velocity set

	// dt(f_i) + c_i Grad(f_i) = -1/tau (f_i - f_eq)
