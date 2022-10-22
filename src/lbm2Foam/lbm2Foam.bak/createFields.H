// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create thermo-hydraulic fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "Reading field U\n" << endl;

	volVectorField U
	(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
	);
	// exract max initial velocity
	scalar UMAX = max(mag(U)).value();
	// create surface scalar field
	surfaceScalarField phi
	(
		IOobject
		(
				"phi",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
		),
		fvc::flux(U)
	);

  Info<< "Reading transportProperties\n" << endl;

  IOdictionary transportProperties
  (
      IOobject
      (
          "transportProperties",
          runTime.constant(),
          mesh,
          IOobject::MUST_READ_IF_MODIFIED,
          IOobject::NO_WRITE
      )
  );

	// kinematic viscosity
	const dimensionedScalar nu
	(
		transportProperties.lookup("nu")
	);

	// ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // Read LBM user options from controlDict

	// read number of directions
	const unsigned int dirNum =
	runTime.controlDict().lookupOrDefault<scalar>("directions", 9);
	// read space discretization parameter
	const scalar deltaX =
	runTime.controlDict().lookupOrDefault<scalar>("deltaX", 1);
	// read time discretization parameter
	scalar deltaT =
	runTime.controlDict().lookupOrDefault<scalar>("deltaT", 1);
	// read flow density
	const scalar rho0 =
	runTime.controlDict().lookupOrDefault<scalar>("density", 1000);
	// read reference pressure
	const scalar pRef =
	runTime.controlDict().lookupOrDefault<scalar>("pressure", 1.01e5);
	// read number of internal pimple loops
	const unsigned int internal =
	runTime.controlDict().lookupOrDefault<scalar>("internal", 1);
	// read system temperature (isothermal model)
	const scalar T0 =
	runTime.controlDict().lookupOrDefault<scalar>("temperature", 900);


	// ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // Dimension sets

	// inverse time dimensions
	const dimensionedScalar IdimTime
	(
		"IdimTime",
		dimensionSet(0,0,-1,0,0,0,0),
		1
	);

	// velocity dimensions
	const dimensionedScalar dimVel
	(
		"dimVel",
		dimensionSet(0,1,-1,0,0,0,0),
		1
	);

	// density dimensions
	const dimensionedScalar dimRho
	(
		"dimRho",
		dimensionSet(1,-3,0,0,0,0,0),
		1
	);

	// pressure dimensions
	const dimensionedScalar dimPres
	(
		"dimPres",
		dimensionSet(1,-1,-2,0,0,0,0),
		1
	);

	// Temperature dimensions
	const dimensionedScalar dimTemp
	(
		"dimTemp",
		dimensionSet(0,0,0,1,0,0,0),
		1
	);

	// ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "Creating LB fields\n" << endl;


	PtrList<volScalarField>				f(dirNum);		// particle distribution
	PtrList<volScalarField>				feq(dirNum);	// equilibrium distribution
	PtrList<volVectorField> 			c(dirNum);		// lattice velocities
	PtrList<surfaceScalarField>		xi(dirNum); 	// lattice velocities surface flux

	// compute velocity set properties
	#include "D2Q9.H"

	// lattice velocity magnitude
	dimensionedScalar cl("cl", dimVel*deltaX/deltaT);
	// inverse of lattice speed of sound squared
	dimensionedScalar ICS2
	(
		"ICS2",
		3/cl/cl
	);
	dimensionedScalar ICS4
	(
		"ICS4",
		ICS2*ICS2
	);

	Info<< "Squared Mach number = " << UMAX*UMAX*ICS2.value() << nl << endl;

	// compute inverse of BGK relaxation time (omega=1/tau)
	scalar tau((0.5 + 3*nu.value()/cl.value()/deltaX)*deltaT);
	dimensionedScalar omega
	(
		"omega",
		IdimTime/tau
	);

	Info<< "BGK (lattice) relaxation time = " << tau/deltaT << nl << endl;

	// particle mass density
	/*
		this is initialize as uniform but computed later
	    according to the initial velocity field
	*/
	volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
				dimRho*rho0
    );
	// incompressible (constant) density field
	volScalarField density("density", rho);
	volScalarField rho_prev("rho_prev", rho);
	// ^ this is needed in simulation initializzation

	// momentum field
	volVectorField momentum(rho*U);

	// pressure field
	/*
		this is related to the density field according to
		the isothermal equation of state ( rho = ICS2*p )
		Actually, since pressure gradients only matters,
		the pressure field can be shifted as:
						rho = rho0 + ICS2*(p - pRef)
	*/
	volScalarField p
  (
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
		dimPres*pRef
  );


	forAll(f, iSet)
	{
		f.set
		(
			iSet,
			new volScalarField
			(
				IOobject
				(
					"f_"+Foam::name(iSet),
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);

		feq.set
		(
			iSet,
			new volScalarField
			(
				IOobject
				(
					"feq_"+Foam::name(iSet),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				0.*dimRho
			)
		);

		c.set
		(
			iSet,
			new volVectorField
			(
				IOobject
				(
					"c_"+Foam::name(iSet),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				cl*lattVect[iSet]
			)
		);

		xi.set
		(
			iSet,
			new surfaceScalarField
			(
				IOobject
				(
					"xi",
					runTime.timeName(),
					mesh,
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE
				),
				fvc::flux(c[iSet]) // this is usually employed for advection terms
			)
		);
	}

// ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	// temperature (constant) field
	dimensionedScalar temperature("T", dimTemp*T0);

	dimensionedScalar dimCv // [= [Cv]/dimTime]
	(
		"dimCv",
		dimensionSet(1,2,-3,-1,0,0,0),
		1
	);

	dimensionedScalar dimHex
	(
		"dimHex",
		dimensionSet(1,2,-3,0,0,0,0),
		1
	);

	dimensionedScalar heatCapacity
	(
		"Cv",
		dimCv*runTime.controlDict().lookupOrDefault<scalar>("heatCapacity", 100)
	);

	dimensionedScalar hexRemoval
	(
		"hex",
		dimHex*runTime.controlDict().lookupOrDefault<scalar>("hex", 7e2)
	);
