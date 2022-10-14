
	// initialize velocity set
	List<scalar> W(dirNum);
	List<vector> lattVect(dirNum);

	W[0] = 4./9.;
	W[1] = 1./9.;
	W[2] = 1./9.;
	W[3] = 1./9.;
	W[4] = 1./9.;
	W[5] = 1./36.;
	W[6] = 1./36.;
	W[7] = 1./36.;
	W[8] = 1./36.;

	lattVect[0] = Foam::vector( 0, 0, 0);
	lattVect[1] = Foam::vector( 1, 0, 0);
	lattVect[2] = Foam::vector( 0, 0, 1);
	lattVect[3] = Foam::vector(-1, 0, 0);
	lattVect[4] = Foam::vector( 0, 0,-1);
	lattVect[5] = Foam::vector( 1, 0, 1);
	lattVect[6] = Foam::vector(-1, 0, 1);
	lattVect[7] = Foam::vector(-1, 0,-1);
	lattVect[8] = Foam::vector( 1, 0,-1);
