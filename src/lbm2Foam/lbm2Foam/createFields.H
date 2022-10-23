


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
	// to compute Mach number
	scalar UMAX = max(mag(U)).value();


  Info<< "Reading transportProperties\n" << endl;

  const IOdictionary transportProperties
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
	dimensionedScalar nu
	(
		transportProperties.lookup("nu")
	);


  // Read LBM user options from controlDict

	// read initialization type
	const bool initScheme =
	runTime.controlDict().lookupOrDefault<bool>("initScheme", false);
	// read space discretization parameter
	const scalar deltaX =
	runTime.controlDict().lookupOrDefault<scalar>("deltaX", 1);
	// read time discretization parameter
	const scalar deltaT =
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



	// Create dimension unities

	// inverse time dimensions
	dimensionedScalar IdimTime
	(
		"IdimTime",
		dimensionSet(0,0,-1,0,0,0,0),
		1
	);

	// velocity dimensions
	dimensionedScalar dimVel
	(
		"dimVel",
		dimensionSet(0,1,-1,0,0,0,0),
		1
	);

	// density dimensions
	dimensionedScalar dimRho
	(
		"dimRho",
		dimensionSet(1,-3,0,0,0,0,0),
		1
	);

	// pressure dimensions
	dimensionedScalar dimPres
	(
		"dimPres",
		dimensionSet(1,-1,-2,0,0,0,0),
		1
	);


	// compute velocity set properties
	#include "readVelocitySet.H"

  Info<< "Creating LBM fields\n" << endl;

	PtrList<volScalarField>				f(dirNum);			// particle distribution
	PtrList<volScalarField>				feq(dirNum);		// equilibrium distribution
	PtrList<volVectorField> 			c(dirNum);			// lattice velocities
	PtrList<volScalarField> 			cDotU(dirNum);  // scalar product between U end c
	PtrList<volScalarField> uEqFactor(dirNum);  	// feq_i = (<-this) * sum_j[ f_j ]
	PtrList<surfaceScalarField>		xi(dirNum); 		// lattice velocities surface flux


	// directions velocity component magnitude
	dimensionedScalar cl(dimVelocity, deltaX/deltaT);
	// save inverse of speed of sound squared
	const dimensionedScalar ICS2
	(
		"ICS2",
		lICS2/cl/cl
	);

	Info<< "Squared Mach number = " << UMAX*UMAX*ICS2.value() << nl << endl;

	// compute inverse of BGK relaxation time (omega=1/tau)
	const scalar tau((0.5 + lICS2*nu.value()/cl.value()/deltaX)*deltaT);
	const dimensionedScalar omega
	(
		"omega",
		IdimTime/tau
	);

	Info<< "BGK (lattice) relaxation time = " << tau/deltaT << nl << endl;

	/* Particle mass density
		This is initialize as uniform but computed later
	  according to the initial velocity field (if required)
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
		dimensionedScalar(dimDensity, rho0)
  );
	// incompressible density field
	volScalarField density("density", rho);
	volScalarField rho_prev("rho_prev", rho);
	// ^ this is needed in simulation initializzation

	/* Pressure field
		This is related to the density field according to
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
		dimensionedScalar(dimPres*pRef)
  );

	// momentum field
	volVectorField momentum("momentum", rho*U);

	// initializzation loop over velocity set
	forAll(f, dI)
	{

		// Read particle distribution
		f.set
		(
			dI,
			new volScalarField
			(
				IOobject
				(
					"f_"+Foam::name(dI),
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);

		// Assign velocity set
		c.set
		(
			dI,
			new volVectorField
			(
				IOobject
				(
					"c_"+Foam::name(dI),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedVector(cl*lattVect[dI])
			)
		);

		/* Save scalar product between lattice velocities
		 	 and velocity field: cDotU = ci.U
		*/
		cDotU.set
		(
			dI,
			new volScalarField
			(
				IOobject
				(
					"cDotU_"+Foam::name(dI),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				(c[dI] & U)
			)
		);

		/* Save initial equilibrium distribution factor:
			The equilibrium distribution in each direction
			is just proportional to the macroscopic density
			field 		feq = uEqFactor * rho ,  where
				uEqFactor =
		wi (1+(ci.U)( 1 + 0.5(ci.U)/cs2)/cs2 - 0.5 U2/cs2 )
		*/
		uEqFactor.set
		(
			dI,
			new volScalarField
			(
				IOobject
				(
					"uEf_"+Foam::name(dI),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				W[dI]*( 1.0
							+ ICS2*cDotU[dI]*(1. + 0.5*ICS2*cDotU[dI])
							- 0.5*ICS2*(U&U))
			)
		);

		// Assign equilibrium distribution values
		feq.set
		(
			dI,
			new volScalarField
			(
				IOobject
				(
					"feq_"+Foam::name(dI),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				uEqFactor[dI]*rho
			)
		);

		// Assign velocity set surface flux values
		xi.set
		(
			dI,
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
				fvc::flux(c[dI])
			)
		);

		// Initialize particle equilibrium distribution
		f[dI] = feq[dI];

	} // end of loop over velocity set
