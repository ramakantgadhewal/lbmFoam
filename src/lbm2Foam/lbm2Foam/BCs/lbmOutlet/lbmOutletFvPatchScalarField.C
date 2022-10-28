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

#include "lbmOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lbmOutletFvPatchScalarField::lbmOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    dI_(Foam::label(0)),
    rhoName_("rho")

{
  #include "computeEquilibriumFactor.H"

  // determine direction
  Field<scalar> relDir(patch().nf()&ci);
  forAll(relDir, fI)
  {
    if(relDir[fI] < 0)
    {
      relDir[fI] = 1.;
    }
    else
    {
      relDir[fI] = 0.;
    }
  }
  orientation_ = relDir;

  // extract outlet pressure
  rhoOut_ = static_cast<Field<scalar>>
  (
    patch().lookupPatchField<volScalarField, scalar>(rhoName_)
  );

}


Foam::lbmOutletFvPatchScalarField::lbmOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    dI_(dict.lookup<Foam::label>("index")),
    rhoName_(dict.lookupOrDefault<word>("moment-0","rho"))
{
  #include "computeEquilibriumFactor.H"

  // determine direction
  Field<scalar> relDir(patch().nf()&ci);
  forAll(relDir, fI)
  {
    if(relDir[fI] < 0)
    {
      relDir[fI] = 1.;
    }
    else
    {
      relDir[fI] = 0.;
    }
  }
  orientation_ = relDir;

  // extract outlet pressure
  rhoOut_ = static_cast<Field<scalar>>
  (
    patch().lookupPatchField<volScalarField, scalar>(rhoName_)
  );

}


Foam::lbmOutletFvPatchScalarField::lbmOutletFvPatchScalarField
(
    const lbmOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    dI_(ptf.dI_),
    rhoName_(ptf.rhoName_),
    orientation_(ptf.orientation_),
    uCEqFactorI_(ptf.uCEqFactorI_),
    rhoOut_(ptf.rhoOut_)
{}


Foam::lbmOutletFvPatchScalarField::lbmOutletFvPatchScalarField
(
    const lbmOutletFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    dI_(tppsf.dI_),
    rhoName_(tppsf.rhoName_),
    orientation_(tppsf.orientation_),
    uCEqFactorI_(tppsf.uCEqFactorI_),
    rhoOut_(tppsf.rhoOut_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::lbmOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    // current component field
    const tmp<Field<scalar>>& fiC = patchInternalField();

    // current component field gradient
    const tmp<Field<vector>>&
    gradfiC =
    patch().lookupPatchField<volVectorField, vector>
    (
      "gradf_"+name(dI_)
    ).patchInternalField();

    // equilibrium distribution
    const tmp<Field<scalar>>&
    fieqC = patch().lookupPatchField<volScalarField, scalar>
    (
      "feq_"+name(dI_)
    ).patchInternalField();

    // extrapolate if unkown distribution
    const tmp<Field<scalar>>&
    dfiC = gradfiC & (patch().Cf() - patch().Cn());
    const Field<scalar>& fiB = static_cast<Field<scalar>>
    (
        fiC
      + orientation_*dfiC
      + (1-orientation_)*(uCEqFactorI_*rhoOut_ - fieqC)
    );
  	operator==
    (
        fiB
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::lbmOutletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

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
        lbmOutletFvPatchScalarField
    );

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
