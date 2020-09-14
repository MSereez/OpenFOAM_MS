/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::skinfrictionCoeff

Group
    grpForcesFunctionObjects

Description
    Computes the skin friction Coeff 


SourceFiles
    skinfrictionCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_skinfrictionCoeff_H
#define functionObjects_skinfrictionCoeff_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class skinfrictionCoeff Declaration
\*---------------------------------------------------------------------------*/

class skinfrictionCoeff
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data
	scalar Uref_;

	scalar rhoRef_;
        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the shear-stress
        void calcskinFriction
        (
            const volSymmTensorField& Reff,
            volVectorField& skinFriction
        );


public:

    //- Runtime type information
    TypeName("skinfrictionCoeff");


    // Constructors

        //- Construct from Time and dictionary
        skinfrictionCoeff
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        skinfrictionCoeff(const skinfrictionCoeff&) = delete;

        //- No copy assignment
        void operator=(const skinfrictionCoeff&) = delete;


    //- Destructor
    virtual ~skinfrictionCoeff() = default;


    // Member Functions

        //- Read the skinfrictionCoeff data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
