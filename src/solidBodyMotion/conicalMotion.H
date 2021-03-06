/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::solidBodyMotionFunctions::conicalMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function.

    The conical motion of an Aerospace body involves rotating around a velocity vector with a given 
    angular velocity.. The solid body motion does not change any properties such as angle of attack or 
    side slip angle.  The appropriate wind tunnel testing equivalent to conical motion is rotary-balance tests.

    Edited by: Mohamed Sereez, demontfortuniversity, Leicester, UK 

SourceFiles
    conicalMotion.C

   

\*---------------------------------------------------------------------------*/

#ifndef conicalMotion_H
#define conicalMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "Function1.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class conicalMotion Declaration
\*---------------------------------------------------------------------------*/

class conicalMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Origin of the axis
        const vector origin_;

        //- Axis vector (velocity vector)
        const vector rotationAxis_;
        //- rotary Angular velocty (rad/sec)
        autoPtr<Function1<scalar>> rotationOmega_;

    // Private Member Functions

        //- Disallow copy construct
        conicalMotion(const conicalMotion&);

        //- Disallow default bitwise assignment
        void operator=(const conicalMotion&);


public:

    //- Runtime type information
    TypeName("conicalMotion");


    // Constructors

        //- Construct from components
        conicalMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new conicalMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~conicalMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
