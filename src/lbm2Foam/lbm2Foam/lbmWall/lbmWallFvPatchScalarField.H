/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::lbmWallFvPatchScalarField

Description
    Boundary condition which switches between maximum and minimum specified
    value if the weigthed average is higher than a threshold specified by the
    user.

SourceFiles
    lbmWallFvPatchScalarField.C

See also
    Foam::fixedValueFvPatchField

\*---------------------------------------------------------------------------*/

#ifndef lbmWallFvPatchScalarField_H
#define lbmWallFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class lbmWallFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class lbmWallFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
  // Private Data

    // wall velocity
    vector uW_;

    // lattice velocity index
    Foam::label dI_;

    // external owner information names
    word rhoName_;
    word feqName_;
    word fName_;

    // velocity set data (written by the default constructor)
    scalar  uwEqFactorI_;

public:

    //- Runtime type information
    TypeName("lbmWall");


    // Constructors

        //- Construct from patch and internal field
        lbmWallFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        lbmWallFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given lbmWallFvPatchScalarField
        //  onto a new patch
        lbmWallFvPatchScalarField
        (
            const lbmWallFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        lbmWallFvPatchScalarField
        (
            const lbmWallFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        lbmWallFvPatchScalarField
        (
            const lbmWallFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new lbmWallFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return wall velocity
            vector wallVelocity() const
            {
                return uW_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
