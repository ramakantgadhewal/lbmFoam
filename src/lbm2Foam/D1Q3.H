	
	// initialize velocity set
	List<scalar> W(3);
	List<vector> lattVect(3);
	
	W[0] = 2./3.;
	W[1] = 1./6.;
	W[2] = 1./6.;
	
	lattVect[0] = Foam::vector( 0, 0, 0);
	lattVect[1] = Foam::vector( 1, 0, 0);
	lattVect[2] = Foam::vector(-1, 0, 0);
	