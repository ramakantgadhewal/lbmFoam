
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

	// Kinematic viscosity
	dimensionedScalar nu
	(
		transportProperties.lookup("nu")
	);

// * * * * * * * * * * * * * * * * * LBM fields * * * * * * * * * * * * * * * /

	// - read initialization type
	const bool initScheme =
	runTime.controlDict().lookupOrDefault<bool>("initScheme", false);
	// - read space discretization parameter
	const scalar deltaX =
	runTime.controlDict().lookupOrDefault<scalar>("deltaX", 1);
	// - read time discretization parameter
	const scalar deltaT =
	runTime.controlDict().lookupOrDefault<scalar>("deltaT", 1);
	// - read flow density
	const scalar rho0 =
	runTime.controlDict().lookupOrDefault<scalar>("density", 1000);
	// - read reference pressure
	const scalar pRef =
	runTime.controlDict().lookupOrDefault<scalar>("pressure", 1.01e5);
	// - read number of internal pimple loops
	const label internal =
	runTime.controlDict().lookupOrDefault<scalar>("internal", 1);
	// BGK relaxation time option
	scalar tau =
	runTime.controlDict().lookupOrDefault<scalar>("BGKtau", 0);



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


	// Compute velocity set properties
	#include "readVelocitySet.H"

  Info<< "Creating LBM fields\n" << endl;

	PtrList<volScalarField>				f(dirNum);			// particle distribution
	PtrList<volScalarField>				feq(dirNum);		// equilibrium distribution
	PtrList<volVectorField> 			c(dirNum);			// lattice velocities
	PtrList<volScalarField> cDotUBycs2(dirNum);  	// scalar product between U end c
	PtrList<volScalarField> uEqFactor(dirNum);  	// feq_i = (<-this) * sum_j[ f_j ]
	PtrList<surfaceScalarField>		xi(dirNum); 		// lattice velocities surface flux
	PtrList<volVectorField> 		gradf(dirNum);  	// needed for inlet&outlet BCs


	// Directions velocity component magnitude
	dimensionedScalar cl(dimVelocity, deltaX/deltaT);

	// Write computed magnitude in velocity set dictionary
	dictionary& cCompDict = velocitySet.subDict("c");
	cCompDict.set("value", cl.value());
	velocitySet.regIOobject::write();

	// Save inverse of speed of sound squared
	const dimensionedScalar ICS2
	(
		"ICS2",
		lICS2/cl/cl
	);

	Info<< "Squared Mach number = " << UMAX*UMAX*ICS2.value() << nl << endl;

	// Compute inverse of BGK relaxation time (omega=1/tau)
	if(tau==0)
	{
		tau = (0.5 + lICS2*nu.value()/cl.value()/deltaX);
	}
	const dimensionedScalar omega
	(
		"omega",
		IdimTime/tau/deltaT
	);

	Info<< "BGK (lattice) relaxation time = " << tau << nl << endl;

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
	// Incompressible density field
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
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
    ),
    mesh
  );

	// Momentum field
	volVectorField momentum("momentum", rho*U);

	// Initializzation loop over velocity set
	forAll(f, dI)
	{
		Info<<"Read particle field component: " << dI << endl;
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
				cl*lattVect[dI]
			)
		);

		/* Save scalar product between lattice velocities
		 	 and velocity field: cDotU = ci.U
		*/
		cDotUBycs2.set
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
				ICS2*(c[dI] & U)
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
							+ cDotUBycs2[dI]*(1. + 0.5*cDotUBycs2[dI])
							- 0.5*ICS2*(U&U) )
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

		gradf.set
		(
			dI,
			new volVectorField
			(
				IOobject
				(
					"gradf_"+Foam::name(dI),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				fvc::grad(f[dI])
			)
		);

	} // end of loop over velocity set

	// separate loop to initialize boundary conditions
	forAll(f, dI)
	{
		f[dI].correctBoundaryConditions();
	}
