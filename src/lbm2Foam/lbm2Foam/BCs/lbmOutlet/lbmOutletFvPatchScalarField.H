/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-----------------------------------------7--------------------------------------
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
    Foam::lbmOutletFvPatchScalarField

Description
    Boundary condition which switches between maximum and minimum specified
    value if the weigthed average is higher than a threshold specified by the
    user.

SourceFiles
    lbmOutletFvPatchScalarField.C

See also
    Foam::fixedValueFvPatchField

\*---------------------------------------------------------------------------*/

#ifndef lbmOutletFvPatchScalarField_H
#define lbmOutletFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class lbmOutletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class lbmOutletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
  // Private Data

    // lattice velocity indices
    Foam::label dI_;            // current

    // density field name
    word rhoName_;

    // tell if distribution component is exiting the wall
    Field<scalar> orientation_;

    // velocity set data (written by the default constructor)
    Field<scalar>  uCEqFactorI_;

    // outlet pressure
    Field<scalar> rhoOut_;


public:

    //- Runtime type information
    TypeName("lbmWall");


    // Constructors

        //- Construct from patch and internal field
        lbmOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        lbmOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given lbmOutletFvPatchScalarField
        //  onto a new patch
        lbmOutletFvPatchScalarField
        (
            const lbmOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        lbmOutletFvPatchScalarField
        (
            const lbmOutletFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        lbmOutletFvPatchScalarField
        (
            const lbmOutletFvPatchScalarField&,
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
                new lbmOutletFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

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
