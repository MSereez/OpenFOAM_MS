/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "forceCoeffs.H"
#include "dictionary.H"
#include "Time.H"
#include "Pstream.H"
#include "IOmanip.H"
#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceCoeffs, 0);
    addToRunTimeSelectionTable(functionObject, forceCoeffs, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forceCoeffs::createFiles()
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !coeffFilePtr_.valid())
    {
        coeffFilePtr_ = createFile("coefficient");
        writeIntegratedHeader("Coefficients", coeffFilePtr_());
    }
}


void Foam::functionObjects::forceCoeffs::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, "Body Force coefficients");
    writeHeaderValue(os, "dragDir", coordSys_.e1());
    writeHeaderValue(os, "sideDir", coordSys_.e2());
    writeHeaderValue(os, "liftDir", coordSys_.e3());
    writeHeaderValue(os, "rollAxis", coordSys_.e1());
    writeHeaderValue(os, "pitchAxis", coordSys_.e2());
    writeHeaderValue(os, "yawAxis", coordSys_.e3());
    writeHeaderValue(os, "magUInf", magUInf_);
    writeHeaderValue(os, "lRef", lRef_);
    writeHeaderValue(os, "Aref", Aref_);
    writeHeaderValue(os, "CofR", coordSys_.origin());
    writeHeaderValue(os, "CofBodyR", bodySys_.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "Cx");
    writeTabbed(os, "Cy");
    writeTabbed(os, "Cz");
    writeTabbed(os, "CMx");
    writeTabbed(os, "CMy");
    writeTabbed(os, "CMz");
    os  << endl;
}


void Foam::functionObjects::forceCoeffs::writeIntegratedData
(
    const word& title,
    const List<Field<scalar>>& coeff
) const
{
    if (!log)
    {
        return;
    }

    const scalar pressure = sum(coeff[0]);
    const scalar viscous = sum(coeff[1]);
    const scalar porous = sum(coeff[2]);
    const scalar total = pressure + viscous + porous;

    Info<< "        " << title << "       : " << total << token::TAB
        << '('
        << "pressure: " << pressure << token::TAB
        << "viscous: " << viscous;

    Info<< ')' << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffs::forceCoeffs
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool readFields
)
:
    forces(name, runTime, dict, false),
    magUInf_(Zero),
    lRef_(Zero),
    lRefRollYaw_(Zero),
    Aref_(Zero),
    coeffFilePtr_()
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict, "liftDir", "dragDir");
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forceCoeffs::read(const dictionary& dict)
{
    forces::read(dict);

    // Free stream velocity magnitude
    dict.readEntry("magUInf", magUInf_);

    // If case is compressible we must read rhoInf (store in rhoRef_) to
    // calculate the reference dynamic pressure
    // Note: for incompressible, rhoRef_ is already initialised
    if (rhoName_ != "rhoInf")
    {
        dict.readEntry("rhoInf", rhoRef_);
    }

    // Reference length and area scales
    dict.readEntry("lRef", lRef_);
    dict.readEntry("lRefRollYaw", lRefRollYaw_);
    dict.readEntry("Aref", Aref_);


    return true;
}


bool Foam::functionObjects::forceCoeffs::execute()
{
    forces::calcForcesMoment();
    forces::getQuaternionBodyAxis();
    forces::calcBodyForcesMoment();

    createFiles();

    scalar pDyn = 0.5*rhoRef_*magUInf_*magUInf_;


	List<Field<scalar>> coeffs(6);
        coeffs[0].setSize(nBin_);
        coeffs[1].setSize(nBin_);
        coeffs[2].setSize(nBin_);
	coeffs[3].setSize(nBin_);
	coeffs[4].setSize(nBin_);
	coeffs[5].setSize(nBin_);
	

	scalar FX = bodyForces_[0];
	scalar FY = bodyForces_[1];
	scalar FZ = bodyForces_[2];
	scalar MX = bodyMoment_[0];
	scalar MY = bodyMoment_[1];
	scalar MZ = bodyMoment_[2];

        coeffs[0] = (FX)/(Aref_*pDyn);
        coeffs[1] = (FY)/(Aref_*pDyn);
        coeffs[2] = (FZ)/(Aref_*pDyn);
        coeffs[3] = (MX)/(Aref_*lRefRollYaw_*pDyn);
        coeffs[4] = (MY)/(Aref_*lRefRollYaw_*pDyn);
        coeffs[5] = (MZ)/(Aref_*lRef_*pDyn);




        scalar CX = sum(coeffs[0]);
        scalar CY = sum(coeffs[1]);
        scalar CZ = sum(coeffs[2]);
        scalar CMX =sum(coeffs[3]);
        scalar CMY =  sum(coeffs[4]) ;
	scalar CMZ = sum(coeffs[5]);
	

        Log << type() << " " << name() << " write:" << nl
            << "    CX    = " << CX << nl
            << "    CY    = " << CY << nl
            << "    CZ    = " << CZ << nl
            << "    CMX = " << CMX << nl
            << "    CMY = " << CMY << nl
	    << "    CMZ = " << CMZ << nl 
	    << "    Body X=" <<  quaternBX_ << nl 
	    << "    Body Y=" <<  quaternBY_ << nl 
	    << "    Body Z=" <<  quaternBZ_ << endl; 


    if (writeToFile())
    {
        writeCurrentTime(coeffFilePtr_());
        coeffFilePtr_()
        << tab << CX << tab  << CY
            << tab << CZ << tab << CMX << tab << CMY << tab << CMZ << endl;
    }



    return true;
}


bool Foam::functionObjects::forceCoeffs::write()
{

    return true;
}


// ************************************************************************* //
