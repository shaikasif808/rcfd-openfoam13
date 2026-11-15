/*---------------------------------------------------------------------------*\
Class
    Foam::functionObjects::speciesReference

Description
    Reference for the rCFD lattice tutorial: the species A, B, C, D, E with
    the reactions of rCFD/rCFD_user_run.h, solved by OpenFOAM's own finite
    volume discretisation on the flux phi of the solver (use the solver
    module "functions" for a frozen flow).

        A -> B      k_AB_wall, only in the cells at the patch wallPatch
        B -> C      k_BC
        C -> D      k_CD
        B + B -> E  k_BE Y_B^2

    Implicit Euler; the sources are linear in the unknown and taken
    implicitly, the species are solved one after the other.

Usage
    speciesReference1
    {
        type        speciesReference;
        libs        ("libspeciesReference.so");
        D           6.578e-7;
        k_AB_wall   1e6;
        k_BC        10;
        k_CD        1;
        k_BE        55.125;
        wallPatch   reactingWall;
        inletPatch  inlet;
        outletPatch outlet;
    }
\*---------------------------------------------------------------------------*/

#ifndef speciesReference_H
#define speciesReference_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "OFstream.H"

namespace Foam
{
namespace functionObjects
{

class speciesReference
:
    public fvMeshFunctionObject
{
    // Private Data

        scalar D_, kAB_, kBC_, kCD_, kBE_;
        word wallPatch_, inletPatch_, outletPatch_, phiName_;

        //- A, B, C, D, E
        PtrList<volScalarField> Y_;

        //- 1 in the cells at the reacting wall
        autoPtr<volScalarField::Internal> wall_;

        autoPtr<OFstream> monitor_;


public:

    TypeName("speciesReference");

    speciesReference(const word& name, const Time& runTime, const dictionary&);

    virtual ~speciesReference();

    virtual bool read(const dictionary&);

    virtual wordList fields() const;

    virtual bool execute();

    virtual bool write();
};

}
}

#endif
