
  // 0. read velocity set data dictionary
  const dictionary& velocitySet =
  db().lookupObject<IOdictionary>
  (
    "velocitySetData"
  );

  // - inverse of squared lattice speed of sound
  const scalar& lICS2
  (
    velocitySet.lookup<scalar>("lICS2")
  );

  // - weights MCD
  const scalar& wMCD
  (
    velocitySet.lookup<scalar>("MCD")
  );

  // - velocity component magnitude
  const scalar& cl
  (
    velocitySet.subDict("c").lookup<scalar>("value")
  );

  // 1. Read velocity patch data
  const tmp<Field<vector>>& UW =
    patch().lookupPatchField<volVectorField, vector>("U");

  const vector& uW = average(static_cast<Field<vector>>(UW));
  // !! ^^ This assumes a single velocity assigned
  // TODO: needs a check

  // 2. Compute lattice parameters

  // - speed of sound
  const scalar& ics2(lICS2/cl/cl);

  // - lattice velocity weigth
  const scalar& wi
  (
    velocitySet.lookup<scalar>
    (
      "wNum_"+Foam::name(dI_)
    )/wMCD
  );

  // - lattice velocity
  const vector& ci
  (
    velocitySet.lookup<vector>
    (
      "lattVect_"+Foam::name(dI_)
    )*cl
  );


  // 3. save wall equilibrium velocity factor
  const scalar& cDotUwBycs2(ics2*(ci & uW));
  uInEqFactorI_ = wi*
                  (
                    1.
                  + cDotUwBycs2*(1+0.5*cDotUwBycs2)
                  - 0.5*ics2*(uW & uW)
                  );
