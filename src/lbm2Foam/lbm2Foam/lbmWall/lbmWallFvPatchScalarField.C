/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "lbmWallFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lbmWallFvPatchScalarField::lbmWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    uW_(Foam::vector(0.,0.,0.)),
    dI_(Foam::label(0)),
    rhoName_("rho"),
    feqName_("feq_0"),
    fName_("f_0")

{
  #include "computeEquilibriumFactor.H"

  // determine direction
  Field<scalar> relDir(patch().nf()&ci);
  forAll(relDir, fI)
  {
    if(relDir[fI] > 0)
    {
      relDir[fI] = 0.;
    }
    else
    {
      relDir[fI] = 1.;
    }
  }
  orientation_ = relDir;

}


Foam::lbmWallFvPatchScalarField::lbmWallFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    uW_(dict.lookupOrDefault<vector>("wallVelocity",Foam::vector(0,0,0))),
    dI_(dict.lookup<Foam::label>("index")),
    rhoName_(dict.lookupOrDefault<word>("moment-0","rho"))
{
  #include "computeEquilibriumFactor.H"

  // determine direction
  Field<scalar> relDir(patch().nf()&ci);
  forAll(relDir, fI)
  {
    if(relDir[fI] > 0)
    {
      relDir[fI] = 0.;
    }
    else
    {
      relDir[fI] = 1.;
    }
  }
  orientation_ = relDir;

  // set field names
  feqName_  = "feq_"+Foam::name(dI_);
  fName_    = "f_"+Foam::name(dI_);

  Info<<"Equilibrium factor = " << uwEqFactorI_
      << " - uW = " << uW_
      << " - index: " << dI_ << endl;
}


Foam::lbmWallFvPatchScalarField::lbmWallFvPatchScalarField
(
    const lbmWallFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    uW_(ptf.uW_),
    dI_(ptf.dI_),
    rhoName_(ptf.rhoName_),
    feqName_(ptf.feqName_),
    fName_(ptf.fName_),
    orientation_(ptf.orientation_),
    uwEqFactorI_(ptf.uwEqFactorI_)
{}


Foam::lbmWallFvPatchScalarField::lbmWallFvPatchScalarField
(
    const lbmWallFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    uW_(tppsf.uW_),
    dI_(tppsf.dI_),
    rhoName_(tppsf.rhoName_),
    feqName_(tppsf.feqName_),
    fName_(tppsf.fName_),
    orientation_(tppsf.orientation_),
    uwEqFactorI_(tppsf.uwEqFactorI_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::lbmWallFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // access the owner internal field
    const Field<scalar>& fiC = patchInternalField();
    Info << "fiC = " << fiC << endl;

    // extract internal field next to the current patch

    // equilibrium distribution
    const Field<scalar>& fieqC =
      patch().lookupPatchField<volScalarField, scalar>
      (
        feqName_
      ).patchInternalField();

    // particle density
    const Field<scalar>& rhoC =
      this->patch().lookupPatchField<volScalarField, scalar>
      (
        rhoName_
      ).patchInternalField();

    Info << "rhoC = " << rhoC <<endl;
    // equilibrium distribution at the boundary
    const Field<scalar>& fieqB = uwEqFactorI_*rhoC;
    Info << "fieqB = " << fieqB << endl;

    // set the boundary value
  	operator==
    (
      //orientation_*( fieqB + (fieqC-fiC) )
      fieqB
    );

    fixedValueFvPatchScalarField::updateCoeffs();
    Info<<"Update coeff completed - index: " << dI_ << endl;
}


void Foam::lbmWallFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeKeyword("wallVelocity")
      << uW_ << token::END_STATEMENT << nl;
    os.writeKeyword("index")
      << dI_ << token::END_STATEMENT << nl;

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        lbmWallFvPatchScalarField
    );

}

// ************************************************************************* //
