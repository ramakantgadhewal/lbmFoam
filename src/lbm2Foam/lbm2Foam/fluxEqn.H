// #######################################################
// LBM for one group neutron diffusion
// #######################################################

	forAll(f_flux, dI)
	{
		Info << "direction " << dI << nl
				 << "xi = "
			 	 << (average(xi[dI])).value() << nl
			   << "f_flux[dI] = "
				 << (average(f_flux[dI])).value() << nl
				 << "omega_n feq_flux = "
				 << (average(fvc::Sp(omega_n,feq_flux[dI]))).value() << nl
				 << "w M00 flux = "
				 << (W[dI]*average(fvc::Sp(M00, flux))).value() << nl
				 << "w delSource = "
				 << (W[dI]*average(delSource)).value() << nl << endl;
		for(unsigned int ic=0; ic<internal; ++ic)
		{
			fvScalarMatrix fFluxEqn
			(
					fvm::ddt(f_flux[dI])
			  + fvm::div(xi[dI], f_flux[dI])
				+ fvm::Sp(omega_n, f_flux[dI])
				==
					fvc::Sp(omega_n, feq_flux[dI])
				+ W[dI] * fvc::Sp(M00, flux)
				+ W[dI] * delSource
			);

			fFluxEqn.relax();
			fvConstraints.constrain(fFluxEqn);
			fFluxEqn.solve();
			fvConstraints.constrain(f_flux[dI]);
			f_flux[dI].correctBoundaryConditions();
			Info << "f_flux"+Foam::name(dI) << " = "
					 << (average(f_flux[dI])).value() << nl << endl;

		} // end internal cycles loop

	} // end loop over velocity set

	// load flux from distribution
	flux *= 0.;
	forAll(f_flux, dI)
	{
		flux += f_flux[dI];
	}
	Info << "Average neutron flux = "
			 << (average(flux)).value() << " n/(m2 s)" << nl <<endl;
	// load equilibrium distribution from flux
	forAll(f_flux, dI)
	{
		feq_flux[dI] = W[dI]*flux;
	}

	// load power density
	powerDens = SP0*flux;
	Info << "Power = "
			 << (fvc::domainIntegrate(powerDens)).value() << " W" <<  nl <<endl;
	// load temperature variation
	temperature += deltaT/heatCapacity*
									(
										fvc::domainIntegrate(powerDens)
										- hexRemoval
									);
	Info << "Temperature = " << temperature.value() << " K" << nl << endl;
	logT = log(temperature/TrefXS);
	// and update Cross Section data
	#include "updateXS.H"

// #######################################################
// Precursors equations
// #######################################################

	delSource *= 0.;
	forAll(prec, pJ)
	{

		fvScalarMatrix precjEqn
		(
				fvm::ddt(prec[pJ])
			+ fvm::div(phi,prec[pJ])
			==
				fvc::Sp(NSF0_0*beta[pJ],flux)
			- fvm::Sp(lam[pJ],prec[pJ])
		);

		precjEqn.relax();
		fvConstraints.constrain(precjEqn);
		precjEqn.solve();
		fvConstraints.constrain(prec[pJ]);
		prec[pJ].correctBoundaryConditions();

		delSource += lam[pJ]*prec[pJ]*cl;

	}
