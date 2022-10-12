


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
	dimensionedScalar nu
	(
		transportProperties.lookup("nu")
	);



    // Read LBM user options from controlDict

	// read number of directions
	int dirNum =
	runTime.controlDict().lookupOrDefault<scalar>("directions", 3);
	// read space discretization parameter
	scalar deltaX =
	runTime.controlDict().lookupOrDefault<scalar>("deltaX", 1);
	// read time discretization parameter
	scalar deltaT =
	runTime.controlDict().lookupOrDefault<scalar>("deltaT", 1);
	// read flow density
	scalar rho0 =
	runTime.controlDict().lookupOrDefault<scalar>("density", 1000);
	// read reference pressure
	scalar pRef =
	runTime.controlDict().lookupOrDefault<scalar>("pressure", 1.01e5);
	// read number of internal pimple loops
	int internal =
	runTime.controlDict().lookupOrDefault<scalar>("internal", 1);



	// Dimension sets

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



    Info<< "Creating LB fields\n" << endl;


	PtrList<volScalarField>			f(dirNum);		// particle distribution
	PtrList<volScalarField>			feq(dirNum);	// equilibrium distribution
	PtrList<volVectorField> 		c(dirNum);		// lattice velocities
	PtrList<surfaceScalarField>		xi(dirNum); 	// lattice velocities surface flux

	// compute velocity set properties
	#include "D1Q3.H"

	// lattice velocity magnitude
	dimensionedScalar cl(dimVel*deltaX/deltaT);
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
		dimensionedScalar(dimRho*rho0)
    );
	// incompressible density field
	volScalarField density("density", rho);
	volScalarField rho_prev("rho_prev", rho);
	// ^ this is needed in simulation initializzation

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
		dimensionedScalar(dimPres*pRef)
    );

	// momentum field
	volVectorField momentum("momentum", rho*U);


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
				dimensionedScalar(dimRho*Foam::scalar(0.))
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
				dimensionedVector(cl*lattVect[iSet])
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
				//fvc::interpolate(c[iSet]) & mesh.Sf() // Eric choice
				fvc::flux(c[iSet]) // this is usually employed for advection terms
			)
		);
	}
