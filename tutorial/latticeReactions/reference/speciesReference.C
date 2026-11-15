#include "speciesReference.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFields.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(speciesReference, 0);
    addToRunTimeSelectionTable(functionObject, speciesReference, dictionary);
}
}


Foam::functionObjects::speciesReference::speciesReference
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);

    const word names[5] = {"A", "B", "C", "D", "E"};

    wordList types
    (
        mesh_.boundary().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    const label inleti = mesh_.boundaryMesh().findIndex(inletPatch_);
    const label walli = mesh_.boundaryMesh().findIndex(wallPatch_);

    if (inleti < 0 || walli < 0)
    {
        FatalErrorInFunction
            << "Patch " << inletPatch_ << " or " << wallPatch_ << " not found"
            << exit(FatalError);
    }

    types[inleti] = fixedValueFvPatchScalarField::typeName;

    Y_.setSize(5);

    forAll(Y_, i)
    {
        IOobject io
        (
            names[i],
            time_.name(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );

        if (io.headerOk())
        {
            Y_.set(i, new volScalarField(io, mesh_));
        }
        else
        {
            Y_.set
            (
                i,
                new volScalarField(io, mesh_, dimensionedScalar(dimless, 0), types)
            );

            Y_[i].boundaryFieldRef()[inleti] == (i == 0 ? 1.0 : 0.0);
        }
    }

    wall_.set
    (
        new volScalarField::Internal
        (
            IOobject("reactingWallCells", time_.name(), mesh_),
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    const labelUList& faceCells = mesh_.boundary()[walli].faceCells();

    forAll(faceCells, i)
    {
        wall_()[faceCells[i]] = 1;
    }

    Info<< type() << " " << name << ": "
        << returnReduce(label(sum(wall_().primitiveField())), sumOp<label>())
        << " cells at " << wallPatch_ << endl;

    if (Pstream::master())
    {
        monitor_.set(new OFstream(time_.globalPath()/"reference_monitor.out"));
        monitor_() << "# time [s], mean A B C D E in the cells at the outlet, "
            << "mass of A B C D E in the domain per density [m3]" << endl;
    }
}


Foam::functionObjects::speciesReference::~speciesReference()
{}


bool Foam::functionObjects::speciesReference::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    D_ = dict.lookup<scalar>("D");
    kAB_ = dict.lookup<scalar>("k_AB_wall");
    kBC_ = dict.lookup<scalar>("k_BC");
    kCD_ = dict.lookup<scalar>("k_CD");
    kBE_ = dict.lookup<scalar>("k_BE");
    wallPatch_ = dict.lookup<word>("wallPatch");
    inletPatch_ = dict.lookup<word>("inletPatch");
    outletPatch_ = dict.lookup<word>("outletPatch");
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    return true;
}


Foam::wordList Foam::functionObjects::speciesReference::fields() const
{
    return wordList(1, phiName_);
}


bool Foam::functionObjects::speciesReference::execute()
{
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    const dimensionedScalar D("D", dimArea/dimTime, D_);
    const dimensionedScalar rT(dimless/dimTime, 1);

    volScalarField& A = Y_[0];
    volScalarField& B = Y_[1];
    volScalarField& C = Y_[2];
    volScalarField& Dsp = Y_[3];
    volScalarField& E = Y_[4];

    const volScalarField::Internal kWall(kAB_*rT*wall_());

    {
        fvScalarMatrix eqn
        (
            fvm::ddt(A) + fvm::div(phi, A, "div(phi,Yi)")
          - fvm::laplacian(D, A, "laplacian(D,Yi)")
          + fvm::Sp(kWall, A)
        );
        eqn.solve("Yi");
    }
    {
        const volScalarField::Internal kB(kBC_*rT + kBE_*rT*B.internalField());

        fvScalarMatrix eqn
        (
            fvm::ddt(B) + fvm::div(phi, B, "div(phi,Yi)")
          - fvm::laplacian(D, B, "laplacian(D,Yi)")
          + fvm::Sp(kB, B)
          - fvm::Su(kWall*A.internalField(), B)
        );
        eqn.solve("Yi");
    }
    {
        fvScalarMatrix eqn
        (
            fvm::ddt(C) + fvm::div(phi, C, "div(phi,Yi)")
          - fvm::laplacian(D, C, "laplacian(D,Yi)")
          + fvm::Sp(kCD_*rT, C)
          - fvm::Su(kBC_*rT*B.internalField(), C)
        );
        eqn.solve("Yi");
    }
    {
        fvScalarMatrix eqn
        (
            fvm::ddt(Dsp) + fvm::div(phi, Dsp, "div(phi,Yi)")
          - fvm::laplacian(D, Dsp, "laplacian(D,Yi)")
          - fvm::Su(kCD_*rT*C.internalField(), Dsp)
        );
        eqn.solve("Yi");
    }
    {
        fvScalarMatrix eqn
        (
            fvm::ddt(E) + fvm::div(phi, E, "div(phi,Yi)")
          - fvm::laplacian(D, E, "laplacian(D,Yi)")
          - fvm::Su(kBE_*rT*sqr(B.internalField()), E)
        );
        eqn.solve("Yi");
    }

    // Monitor: the same quantities as rCFD_lattice_monitor.out
    const label outi = mesh_.boundaryMesh().findIndex(outletPatch_);
    const labelUList& outCells = mesh_.boundary()[outi].faceCells();
    const scalarField& V = mesh_.V();

    scalar outV = 0;
    FixedList<scalar, 5> out(0.0), tot(0.0);

    forAll(outCells, i)
    {
        outV += V[outCells[i]];

        forAll(Y_, k)
        {
            out[k] += V[outCells[i]]*Y_[k][outCells[i]];
        }
    }

    forAll(Y_, k)
    {
        tot[k] = gSum(V*Y_[k].primitiveField());
        reduce(out[k], sumOp<scalar>());
    }

    reduce(outV, sumOp<scalar>());

    if (Pstream::master())
    {
        monitor_() << time_.value();
        forAll(Y_, k) monitor_() << ' ' << out[k]/outV;
        forAll(Y_, k) monitor_() << ' ' << tot[k];
        monitor_() << endl;
    }

    return true;
}


bool Foam::functionObjects::speciesReference::write()
{
    return true;
}
