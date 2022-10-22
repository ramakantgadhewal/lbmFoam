	// number of precursors groups
	unsigned int NPREC(8);

	// Dimension sets

	// neutron flux dimensions
	const dimensionedScalar dimFlux
	(
		"dimFlux",
		dimensionSet(0,-2,-1,0,0,0,0),
		1
	);

	// precursor concentration dimension
	const dimensionedScalar dimPrec
	(
		"dimPrec",
		dimensionSet(0,-3,0,0,0,0,0),
		1
	);


	Info<< "Reading nuclear properties.\n" << endl;

	IOdictionary nuclearProperties
	(
		IOobject
		(
				"nuclearProperties",
				runTime.constant(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
		)
	);

	// neutron XS data
	const dimensionedScalar // Inverse neutron speed [s/m]
		IV0(nuclearProperties.lookup("IV0"));
	const dimensionedScalar // Diffusion coefficient [1/m]
		D0_0(nuclearProperties.lookup("D0_0"));
	const dimensionedScalar // Absorption XS [1/m]
		A0_0(nuclearProperties.lookup("A0_0"));
	const dimensionedScalar // nu*XSf: (mean neutron gen/fiss)*(fission XS) [1/m]
		NSF0_0(nuclearProperties.lookup("NSF0_0"));
	const dimensionedScalar // power XS [W m2/s]
		SP0_0(nuclearProperties.lookup("SP0_0"));
	// reactivity log coefficients [same dimensions]
	const dimensionedScalar
		alfa_D0(nuclearProperties.lookup("alfa_D0"));
	const dimensionedScalar
		alfa_A0(nuclearProperties.lookup("alfa_A0"));
	const dimensionedScalar
		alfa_NSF0(nuclearProperties.lookup("alfa_NSF0"));
	const dimensionedScalar
		alfa_SP0(nuclearProperties.lookup("alfa_SP0"));
	// XS data eference temperature
	const dimensionedScalar
		TrefXS(nuclearProperties.lookup("TrefXS"));

	Info<< "Creating XS fields.\n" << endl;

	dimensionedScalar logT("logT", log(temperature/TrefXS));
	dimensionedScalar NSF0("NSF0", NSF0_0 + alfa_NSF0*logT);
	dimensionedScalar D0("D0", D0_0 + alfa_D0*logT);
	dimensionedScalar SP0("SP0", SP0_0 + alfa_SP0*logT);

  Info<< "Creating LB neutron diffusion fields.\n" << endl;

	PtrList<volScalarField>			f_flux(dirNum);		// particle distribution
	PtrList<volScalarField>			feq_flux(dirNum);	// equilibrium distribution


	// compute inverse of BGK relaxation time (omega=1/tau)
	scalar tau_n((0.5 + 3*(D0_0/deltaX).value())*deltaT);
	dimensionedScalar omega_n
	(
		"omega_n",
		IdimTime/tau_n
	);

	Info<< "BGK (lattice) neutron relaxation time = " << tau_n/deltaT << nl << endl;

	// neutron density
	/*
		this is initialize as uniform but computed later
	    according to the initial velocity field
	*/
	volScalarField flux
  (
    IOobject
    (
        "flux",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
	);

	forAll(f_flux, iSet)
	{
		f_flux.set
		(
			iSet,
			new volScalarField
			(
				IOobject
				(
					"f_flux_"+Foam::name(iSet),
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);

		feq_flux.set
		(
			iSet,
			new volScalarField
			(
				IOobject
				(
					"feq_flux_"+Foam::name(iSet),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar(dimFlux*Foam::scalar(0.))
			)
		);

	}

	Info<< "Reading precursors data.\n" << endl;

	PtrList<dimensionedScalar> 	lam(NPREC);
	PtrList<dimensionedScalar> 	beta(NPREC);
	dimensionedScalar						betaEff(0.);
	PtrList<volScalarField> 		prec(NPREC);

	IOdictionary delayedNeutronConstants
	(
		IOobject
		(
				"delayedNeutronConstants",
				runTime.constant(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
		)
	);

	forAll(prec, pI)
	{

		// read precursors decay constant
		lam.set
		(
			pI,
			new dimensionedScalar
			(
				delayedNeutronConstants.lookup("lam"+Foam::name(pI+1))
			)
		);

		// read precursors fraction
		beta.set
		(
			pI,
			new dimensionedScalar
			(
				delayedNeutronConstants.lookup("beta"+Foam::name(pI+1))
			)
		);
		betaEff += beta[pI];

		// create precursors concentration fields
		prec.set
		(
			pI,
			new volScalarField
			(
				IOobject
				(
					"prec"+Foam::name(pI+1),
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);

	}

	dimensionedScalar M00("M00", cl*((1-betaEff)*NSF0-(A0_0+alfa_A0*logT)));

	// delayed neutron source
	volScalarField delSource
	(
		"delSource",
		0.*IdimTime*flux
	);

	// power density
	volScalarField powerDens
	(
		"powerDens",
		SP0*flux
	);
