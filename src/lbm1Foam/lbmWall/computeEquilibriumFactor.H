
  // read velocity set data dictionary
  const dictionary& velocitySet =
  db().lookupObject<IOdictionary>
  (
    "velocitySetData"
  );

  const dictionary& controlDict =
    velocitySet.parent().parent();

  // 0. Read discretization parameters

  // read data dictionary

/*
  const dictionary& controlDict =
  db().lookupObject<IOdictionary>
  (
    "controlDict"
  );
*/

  // read space discretization parameter
  const scalar deltaX
  (
    controlDict.lookupOrDefault<scalar>("deltaX", 1)
  );

  // read time discretization parameter
  const scalar deltaT
  (
    controlDict.lookupOrDefault<scalar>("deltaT", 1)
  );


  // 1. Read velocity set data

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


  // 2. Compute lattice parameters

  // speed of sound
  const scalar ics2(lICS2*deltaT*deltaT/deltaX/deltaX);

  // lattice velocity weigth
  const scalar wi
  (
    velocitySet.lookup<scalar>
    (
      "wNum_"+Foam::name(dI_)
    )/wMCD
  );

  // lattice velocity
  const vector ci
  (
    velocitySet.lookup<vector>
    (
      "lattVect_"+Foam::name(dI_)
    )*deltaX/deltaT
  );


  // 3. evaluate equilibrium wall velocity factor
  const scalar xi(ics2*(ci & uW_));
  uwEqFactorI_ = wi*
                  (
                    1.
                  + xi*(1+0.5*xi)
                  - 0.5*ics2*(uW_ & uW_)
                  );
