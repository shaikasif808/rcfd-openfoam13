/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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

#include "BeetstraLimited.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(BeetstraLimited, 0);
    addToRunTimeSelectionTable(dragModel, BeetstraLimited, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::BeetstraLimited::BeetstraLimited
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::BeetstraLimited::~BeetstraLimited()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::BeetstraLimited::CdRe() const
{
    const volScalarField alpha1
    (
        max(interface_.dispersed(), interface_.dispersed().residualAlpha())
    );

    const volScalarField alpha2
    (
        max(interface_.continuous(), interface_.continuous().residualAlpha())
    );

    // OpenFOAM-13 takes Res = max(alpha2*Re, 0) and then evaluates
    // pow(Res, -0.343): a division by zero wherever the two phases have the
    // same velocity (initial field, walls). The limit below removes it;
    // at Res = 1e-3 the drag is the Stokes limit F0 to within 1e-3.
    const volScalarField Res(max(alpha2*interface_.Re(), scalar(1e-3)));

    const volScalarField F0
    (
        "F0",
        10*alpha1/sqr(alpha2) + sqr(alpha2)*(1 + 1.5*sqrt(alpha1))
    );

    const volScalarField F1
    (
        "F1",
        0.413*Res/(24*sqr(alpha2))*(1.0/alpha2
        + 3*alpha1*alpha2 + 8.4*pow(Res, -0.343))
        /(1 + pow(10.0, 3*alpha1)*pow(Res, -(1 + 4*alpha1)/2.0))
    );

    return 24*alpha2*(F0 + F1);
}


// ************************************************************************* //
