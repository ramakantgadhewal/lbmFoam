
  // 0. read velocity set data dictionary
  const dictionary& velocitySet =
  db().lookupObject<IOdictionary>
  (
    "velocitySetData"
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

  // velocity component magnitude
  const scalar cl
  (
    velocitySet.subDict("c").lookup<scalar>("value")
  );


  // 2. Compute lattice parameters

  // speed of sound
  const scalar ics2(lICS2/cl/cl);

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
    )*cl
  );


  // 3. evaluate equilibrium wall velocity factor
  const scalar xi(ics2*(ci & uW_));
  uwEqFactorI_ = wi*
                  (
                    1.
                  + xi*(1+0.5*xi)
                  - 0.5*ics2*(uW_ & uW_)
                  );
