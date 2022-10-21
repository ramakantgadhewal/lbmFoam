/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 23 "/home/mala/OpenFOAM/mala-9/platforms/linux64GccDPInt32Opt/src/solvers/lbmFoam/test/lbm1Foam/periodicLBM/0/U/#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_89843ee9da341cce18cbe91b818ce296e17e42bd
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 40 "/home/mala/OpenFOAM/mala-9/platforms/linux64GccDPInt32Opt/src/solvers/lbmFoam/test/lbm1Foam/periodicLBM/0/U/#codeStream"
// profile parameters declaration
        scalar IwidthU2(4);
        scalar maxU(1);
        scalar xLoc;
        scalar uX;

        const IOdictionary& d = static_cast<const IOdictionary&>
        (
          dict.parent().parent()
        );

        const fvMesh& mesh = refCast<const fvMesh>(d.db());
        // below for patch fields
        //const label id = mesh.boundary().findPatchID("rightWall");
        //const fvPatch& patch = mesh.boundary()[id];

        // initialize field
        vectorField U = volVectorField(mesh.nCells(), vector(0.,0.,0.))

        forAll(mesh.C(), cI)
        {
          xLoc = mesh.C()[cI].x();
          uX = maxU*exp(-xLoc*xLoc*IwidthU2);
          U[cI] = vector(uX,0.,0.);
        }

        writeEntry(os,"",fld);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

