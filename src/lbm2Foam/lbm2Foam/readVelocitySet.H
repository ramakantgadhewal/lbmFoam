/* Read the velocity set dictionary

  This header produces lists of scalars
    - W       : lattice weights
    - lattVect: lattice vectors
  and lattice parameters:
    - dirNum  : (int) number of velocity directions
    - lICS2   : (scalar) inverse of squared speed of sound
  according to the data contained in the dictionary
    constant/velocitySetData
*/

  Info<< "Reading velocity set.\n" << endl;

	IOdictionary velocitySet
	(
		IOobject
		(
			"velocitySetData",
			runTime.constant(),
			mesh,
			IOobject::MUST_READ_IF_MODIFIED,
			IOobject::AUTO_WRITE
		)
	);

  // directions number
	const unsigned int dirNum
	(
		velocitySet.lookup<int>("dirNum")
	);

  // inverse of squared lattice speed of sound
	const scalar lICS2
	(
		velocitySet.lookup<scalar>("lICS2")
	);

  // weights MCD
	const scalar wMCD
	(
		velocitySet.lookup<scalar>("MCD")
	);

  // initialize lists
	List<scalar> W(dirNum);
	List<vector> lattVect(dirNum);

	forAll(W, i)
	{
    // lattice weights
		W[i] = velocitySet.lookup<scalar>
		(
			"wNum_"+Foam::name(i)
		)/wMCD;

    // lattice vectors
		lattVect[i] = velocitySet.lookup<vector>
		(
			"lattVect_"+Foam::name(i)
		);
	}

	Info<< " Number of directions = "
		  << dirNum << nl << endl;
