/*---------------------------------------------------------------------------*\
  rCFD preparation mode inside an OpenFOAM run (OpenFOAM-13)       TUHH 11/2026
\*---------------------------------------------------------------------------*/

#include "rCFDPrep.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmDdt.H"
#include "DynamicList.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvmSup.H"
#include "convectionScheme.H"
#include "fvMatrix.H"
#include "zeroGradientFvPatchFields.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "syncTools.H"
#include "labelIOList.H"
#include "IStringStream.H"
#include "cpuTime.H"

#include "rcfd_of_bridge.h"

#include <fstream>
#include <vector>
#include <utility>
#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(rCFDPrep, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        rCFDPrep,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::rCFDPrep::cellAndCoronaValues
(
    const scalarField& cellValues,
    List<double>& values,
    const label stride,
    const label offset
) const
{
    const label nInt = mesh_.nCells();

    scalarList nbrValues;
    syncTools::swapBoundaryCellList(mesh_, cellValues, nbrValues);

    forAll(cellValues, celli)
    {
        values[stride*celli + offset] = cellValues[celli];
    }

    forAll(coronaFace_, i)
    {
        values[stride*(nInt + i) + offset] = nbrValues[coronaFace_[i]];
    }
}


void Foam::functionObjects::rCFDPrep::open()
{
    if (mesh_.dynamic())
    {
        FatalErrorInFunction
            << "Moving or changing meshes are not supported"
            << exit(FatalError);
    }

    const label nInt = mesh_.nCells();
    const label nIntFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Cell index in the undecomposed mesh
    labelList gid(identityMap(nInt));

    if (Pstream::parRun())
    {
        // A mesh made in parallel (snappyHexMesh -parallel) has no
        // addressing until reconstructPar has been run once
        if
        (
           !IOobject
            (
                "cellProcAddressing",
                mesh_.facesInstance(),
                polyMesh::meshSubDir,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ).headerOk()
        )
        {
            FatalErrorInFunction
                << "No cellProcAddressing in the processor directories "
                << "(instance " << mesh_.facesInstance() << ")." << nl
                << "The mesh was made in parallel: run reconstructPar once "
                << "before foamRun, it reconstructs the mesh and writes the "
                << "addressing that numbers the cells of the database."
                << exit(FatalError);
        }

        gid =
            labelIOList
            (
                IOobject
                (
                    "cellProcAddressing",
                    mesh_.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

        if (gid.size() != nInt)
        {
            FatalErrorInFunction
                << "cellProcAddressing has " << gid.size()
                << " entries, the mesh of this processor has " << nInt
                << " cells" << exit(FatalError);
        }
    }

    const label nGlobal = returnReduce(nInt, sumOp<label>());

    // Values of the cells on the other side of the boundary faces
    labelList nbrGid;
    syncTools::swapBoundaryCellList(mesh_, gid, nbrGid);

    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.cellVolumes();
    const vectorField& Sf = mesh_.faceAreas();

    scalarList nbrV;
    syncTools::swapBoundaryCellList(mesh_, V, nbrV);

    FixedList<scalarList, 3> nbrC;
    for (direction d = 0; d < 3; d++)
    {
        const scalarField Cd(C.component(d));
        syncTools::swapBoundaryCellList(mesh_, Cd, nbrC[d]);
    }

    // Partition faces: the faces of the processor patches
    DynamicList<label> pfFace;
    DynamicList<label> pfHost;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<processorCyclicPolyPatch>(pp))
        {
            FatalErrorInFunction
                << "Patch " << pp.name() << ": cyclic patches are not "
                << "supported" << exit(FatalError);
        }
        else if (isA<processorPolyPatch>(pp))
        {
            const label host =
                refCast<const processorPolyPatch>(pp).neighbProcNo();

            forAll(pp, i)
            {
                pfFace.append(pp.start() + i);
                pfHost.append(host);
            }
        }
        else if (pp.coupled())
        {
            FatalErrorInFunction
                << "Patch " << pp.name() << " of type " << pp.type()
                << ": coupled patches other than processor patches are not "
                << "supported" << exit(FatalError);
        }
    }

    // Corona cells: the distinct cells behind the partition faces, ordered
    // by host and then by global cell index
    std::vector<std::pair<long long, label>> keys(pfFace.size());

    forAll(pfFace, i)
    {
        keys[i].first =
            static_cast<long long>(pfHost[i])*static_cast<long long>(nGlobal)
          + static_cast<long long>(nbrGid[pfFace[i] - nIntFaces]);
        keys[i].second = i;
    }

    std::sort(keys.begin(), keys.end());

    labelList pfCorona(pfFace.size(), -1);
    DynamicList<label> coronaFace;

    for (std::size_t k = 0; k < keys.size(); k++)
    {
        if (k == 0 || keys[k].first != keys[k - 1].first)
        {
            coronaFace.append(pfFace[keys[k].second] - nIntFaces);
        }

        pfCorona[keys[k].second] = coronaFace.size() - 1;
    }

    coronaFace_ = coronaFace;

    const label nExt = coronaFace_.size();
    const label nCells = nInt + nExt;
    const label nFaces = nIntFaces + pfFace.size();

    // Arrays for the library
    std::vector<double> x(3*nCells), vol(nCells), area(3*nFaces);
    std::vector<int> host(nCells), gidAll(nCells), c0(nFaces), c1(nFaces);

    for (label celli = 0; celli < nInt; celli++)
    {
        for (direction d = 0; d < 3; d++)
        {
            x[3*celli + d] = C[celli][d];
        }

        vol[celli] = V[celli];
        host[celli] = Pstream::myProcNo();
        gidAll[celli] = gid[celli];
    }

    forAll(pfFace, i)
    {
        const label bFacei = pfFace[i] - nIntFaces;
        const label cellj = nInt + pfCorona[i];

        // Every face of a corona cell carries the same values
        for (direction d = 0; d < 3; d++)
        {
            x[3*cellj + d] = nbrC[d][bFacei];
        }

        vol[cellj] = nbrV[bFacei];
        host[cellj] = pfHost[i];
        gidAll[cellj] = nbrGid[bFacei];
    }

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    for (label facei = 0; facei < nIntFaces; facei++)
    {
        c0[facei] = owner[facei];
        c1[facei] = neighbour[facei];

        for (direction d = 0; d < 3; d++)
        {
            area[3*facei + d] = Sf[facei][d];
        }
    }

    forAll(pfFace, i)
    {
        const label facei = pfFace[i];
        const label k = nIntFaces + i;

        c0[k] = mesh_.faceOwner()[facei];
        c1[k] = nInt + pfCorona[i];

        // Area vector of a boundary face points out of the own cell
        for (direction d = 0; d < 3; d++)
        {
            area[3*k + d] = Sf[facei][d];
        }
    }

    rcfd_of_mesh_t m;

    m.n_int = nInt;
    m.n_ext = nExt;
    m.n_fint = nIntFaces;
    m.n_fext = pfFace.size();
    m.n_cells_global = nGlobal;
    m.x = x.data();
    m.vol = vol.data();
    m.host = host.data();
    m.gid = gidAll.data();
    m.c0 = c0.data();
    m.c1 = c1.data();
    m.area = area.data();

    rcfd_of_prep_settings_t s;

    rcfd_of_prep_default_settings(&s);

    string phaseList;
    forAll(phaseNames_, phasei)
    {
        phaseList += (phasei ? " " : "") + phaseNames_[phasei];
    }

    s.user_library = userLibrary_.c_str();
    s.n_phases = phaseNames_.size();
    s.phase_names = phaseList.c_str();
    s.norm_phase = normPhase_;

    Info<< type() << " " << name() << ": opening the rCFD library ("
        << nGlobal << " cells, phases: " << phaseList.c_str()
        << ", user library " << userLibrary_.c_str() << ")" << endl;

    rcfd_of_prep_info_t info;

    rcfd_of_prep_open(&m, &s, &info);

    // Boundary patches for the run mode (the counterpart of a face loop over
    // a boundary thread in Fluent): data/topo/patches_<rank>, one line per
    // face with the cell next to it, the face area and the face centre.
    // Every rank lists every patch, with no faces where it has none.
    {
        const fileName patchFile
        (
            "data/topo/patches_" + Foam::name(Pstream::myProcNo())
        );

        std::ofstream os(patchFile.c_str());

        if (!os.good())
        {
            FatalErrorInFunction
                << "Cannot write " << patchFile << exit(FatalError);
        }

        os.precision(17);

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const scalarField magSf(mag(mesh_.faceAreas()));
        const vectorField& Cf = mesh_.faceCentres();

        label nPatches = 0;
        forAll(pbm, patchi)
        {
            if (!pbm[patchi].coupled()) nPatches++;
        }

        os << nPatches << '\n';

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (pp.coupled()) continue;

            const labelUList& faceCells = pp.faceCells();

            os << pp.name() << ' ' << pp.type() << ' ' << pp.size() << '\n';

            forAll(pp, i)
            {
                const label facei = pp.start() + i;

                os << faceCells[i] << ' ' << magSf[facei] << ' '
                   << Cf[facei].x() << ' ' << Cf[facei].y() << ' '
                   << Cf[facei].z() << '\n';
            }
        }
    }

    frameInterval_ = info.frame_interval;
    nFrames_ = info.n_frames;

    Info<< type() << " " << name() << ": " << nFrames_ << " frames of "
        << frameInterval_ << " s from t = " << recordStart_ << " s"
        << " (rCFD_user_prep.h)" << endl;

    opened_ = true;
}


void Foam::functionObjects::rCFDPrep::resetLabels()
{
    forAll(labels_, i)
    {
        volScalarField& X = labels_[i];

        X.primitiveFieldRef() = mesh_.C().primitiveField().component(i % 3);
        X.correctBoundaryConditions();
    }
}


void Foam::functionObjects::rCFDPrep::solveLabels(const label phasei)
{
    const surfaceScalarField& flux =
        mesh_.lookupObject<surfaceScalarField>(fluxNames_[phasei]);

    const bool massFlux = (flux.dimensions() == dimMass/dimTime);

    if (!massFlux && flux.dimensions() != dimVolume/dimTime)
    {
        FatalErrorInFunction
            << "Incompatible dimensions for " << flux.name() << ": "
            << flux.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
    }

    if (massFlux && rhoNames_[phasei] == "none")
    {
        FatalErrorInFunction
            << flux.name() << " is a mass flux: the density field of phase "
            << phaseNames_[phasei] << " is needed (fields { "
            << phaseNames_[phasei] << " { rho ...; } })" << exit(FatalError);
    }

    const bool hasAlpha = (alphaNames_[phasei] != "none");

    const volScalarField* alphaPtr =
        hasAlpha
      ? &mesh_.lookupObject<volScalarField>(alphaNames_[phasei])
      : nullptr;

    const volScalarField* rhoPtr =
        massFlux
      ? &mesh_.lookupObject<volScalarField>(rhoNames_[phasei])
      : nullptr;

    // Continuity error of the time step: taken out of the label equation so
    // that a uniform label stays uniform
    tmp<volScalarField> tcontErr;

    if (hasAlpha && massFlux)
    {
        tcontErr = fvc::ddt(*alphaPtr, *rhoPtr) + fvc::div(flux);
    }
    else if (hasAlpha)
    {
        tcontErr = fvc::ddt(*alphaPtr) + fvc::div(flux);
    }
    else if (massFlux)
    {
        tcontErr = fvc::ddt(*rhoPtr) + fvc::div(flux);
    }
    else
    {
        tcontErr = fvc::div(flux);
    }

    const volScalarField& contErr = tcontErr();

    const dimensionedScalar residualAlpha(dimless, residualAlpha_);

    const Vector<label>& solutionD = mesh_.solutionD();

    // Courant number of the label transport in the cells that hold the phase.
    // The solver may run at a large Courant number (implicit, sub-cycled
    // phase fractions); the labels are then advanced in sub-steps with the
    // flux of the time step, so that their error stays that of a small one.
    typedef DimensionedField<scalar, volMesh> dField;

    const scalar deltaT = time_.deltaTValue();

    tmp<dField> tmNew, tmOld;

    if (hasAlpha && massFlux)
    {
        tmNew = alphaPtr->internalField()*rhoPtr->internalField();
        tmOld =
            alphaPtr->oldTime().internalField()
           *rhoPtr->oldTime().internalField();
    }
    else if (hasAlpha)
    {
        tmNew = new dField("mNew", alphaPtr->internalField());
        tmOld = new dField("mOld", alphaPtr->oldTime().internalField());
    }
    else if (massFlux)
    {
        tmNew = new dField("mNew", rhoPtr->internalField());
        tmOld = new dField("mOld", rhoPtr->oldTime().internalField());
    }
    else
    {
        const IOobject io("mOne", time_.name(), mesh_);

        tmNew = new dField(io, mesh_, dimensionedScalar(dimless, 1));
        tmOld = new dField(io, mesh_, dimensionedScalar(dimless, 1));
    }

    const dField& mNew = tmNew();
    const dField& mOld = tmOld();

    label nSub = 1;
    {
        // Courant number of the velocity that carries the labels. A face flux
        // is divided by the larger content of its two cells (the upwind one
        // at a front): dividing by the content of the cell itself would count
        // the filling of a nearly empty cell as a large Courant number.
        const scalarField& V = mesh_.V();
        const scalarField& m = mNew;
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const scalarField& fluxI = flux.primitiveField();

        scalarField sumFlux(mesh_.nCells(), 0);

        forAll(nei, facei)
        {
            const scalar w =
                mag(fluxI[facei])/max(max(m[own[facei]], m[nei[facei]]), vSmall);

            sumFlux[own[facei]] += w;
            sumFlux[nei[facei]] += w;
        }

        scalarList nbrM;
        syncTools::swapBoundaryCellList(mesh_, m, nbrM);

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];
            const labelUList& faceCells = pp.faceCells();
            const scalarField& fluxP = flux.boundaryField()[patchi];

            forAll(fluxP, i)
            {
                const label celli = faceCells[i];

                const scalar mf =
                    pp.coupled()
                  ? max(m[celli], nbrM[pp.start() + i - mesh_.nInternalFaces()])
                  : m[celli];

                sumFlux[celli] += mag(fluxP[i])/max(mf, vSmall);
            }
        }

        // Where the phase is nearly absent the ratio of flux and content is
        // noise: those cells do not decide the number of sub-steps
        const scalar alphaMin = max(1e-3, residualAlpha_);

        scalar co = 0;

        forAll(V, celli)
        {
            if (hasAlpha && (*alphaPtr)[celli] < alphaMin) continue;

            co = max(co, 0.5*deltaT*sumFlux[celli]/V[celli]);
        }

        reduce(co, maxOp<scalar>());

        if (co > labelMaxCo_)
        {
            nSub = min
            (
                maxLabelSubCycles_,
                label(ceil(co/labelMaxCo_ - 1e-9))
            );
        }

        coFrame_ = max(coFrame_, co);
        nSubFrame_ = max(nSubFrame_, nSub);
    }

    if (nSub > 1)
    {
        // Sub-step k:  (m1 X1 - m0 X0)/dt + div(flux X1) - contErr X1 = 0
        // with the content m = alpha rho interpolated between the old and
        // the new time; (m1 - m0)/dt + div(flux) = contErr in every sub-step
        const dimensionedScalar rDtSub(dimless/dimTime, nSub/deltaT);

        tmp<dField> tmRes;

        if (hasAlpha && massFlux)
        {
            tmRes = residualAlpha_*rhoPtr->internalField();
        }
        else if (hasAlpha)
        {
            tmRes =
                new dField
                (
                    IOobject("mRes", time_.name(), mesh_),
                    mesh_,
                    residualAlpha
                );
        }
        else
        {
            tmRes = 0.0*mNew;
        }

        const dField& mRes = tmRes();

        for (label subi = 0; subi < nSub; subi++)
        {
            const scalar f0 = scalar(subi)/nSub;
            const scalar f1 = scalar(subi + 1)/nSub;

            const dField c0("c0", (mOld + f0*(mNew - mOld) + mRes)*rDtSub);
            const dField c1("c1", (mOld + f1*(mNew - mOld) + mRes)*rDtSub);

            for (direction d = 0; d < 3; d++)
            {
                if (solutionD[d] == -1) continue;

                volScalarField& X = labels_[3*phasei + d];

                IStringStream scheme(divScheme_);

                tmp<fvScalarMatrix> tXEqn
                (
                    fv::convectionScheme<scalar>::New(mesh_, flux, scheme)()
                   .fvmDiv(flux, X)
                  + fvm::Sp(c1 - contErr.internalField(), X)
                  - fvm::Su(c0*X.internalField(), X)
                );

                tXEqn.ref().solve(solverDict_);
            }
        }

        return;
    }

    for (direction d = 0; d < 3; d++)
    {
        if (solutionD[d] == -1) continue;

        volScalarField& X = labels_[3*phasei + d];

        IStringStream scheme(divScheme_);

        tmp<fvScalarMatrix> tXEqn
        (
            fv::convectionScheme<scalar>::New(mesh_, flux, scheme)()
           .fvmDiv(flux, X)
          - fvm::Sp(contErr, X)
        );
        fvScalarMatrix& XEqn = tXEqn.ref();

        if (hasAlpha && massFlux)
        {
            const volScalarField& alpha = *alphaPtr;
            const volScalarField& rho = *rhoPtr;

            XEqn +=
                fvm::ddt(alpha, rho, X)
              + fvm::ddt(residualAlpha_*rho, X)
              - fvc::ddt(residualAlpha_*rho, X);
        }
        else if (hasAlpha)
        {
            XEqn +=
                fvm::ddt(*alphaPtr, X)
              + fvm::ddt(residualAlpha, X)
              - fvc::ddt(residualAlpha, X);
        }
        else if (massFlux)
        {
            XEqn += fvm::ddt(*rhoPtr, X);
        }
        else
        {
            XEqn += fvm::ddt(X);
        }

        XEqn.solve(solverDict_);
    }
}


void Foam::functionObjects::rCFDPrep::handOverFrame()
{
    const label nCells = mesh_.nCells() + coronaFace_.size();
    const label nPhases = phaseNames_.size();

    // vof[cell*nPhases + phase], origin[(phase*nCells + cell)*3 + direction]
    List<double> vof(nCells*nPhases), origin(nPhases*3*nCells);

    const boundBox& bb = mesh_.bounds();

    forAll(phaseNames_, phasei)
    {
        if (alphaNames_[phasei] == "none")
        {
            for (label i = 0; i < nCells; i++)
            {
                vof[i*nPhases + phasei] = 1;
            }
        }
        else
        {
            const volScalarField& alpha =
                mesh_.lookupObject<volScalarField>(alphaNames_[phasei]);

            cellAndCoronaValues(alpha.primitiveField(), vof, nPhases, phasei);

            for (label i = 0; i < nCells; i++)
            {
                double& a = vof[i*nPhases + phasei];
                a = min(max(a, 0.0), 1.0);
            }
        }

        const label offset = phasei*3*nCells;

        for (direction d = 0; d < 3; d++)
        {
            cellAndCoronaValues
            (
                labels_[3*phasei + d].primitiveField(),
                origin,
                3,
                offset + d
            );

            // An origin cannot lie outside the domain (a limited scheme may
            // overshoot by a small amount next to a wall)
            for (label i = 0; i < nCells; i++)
            {
                double& x = origin[offset + 3*i + d];
                x = min(max(x, bb.min()[d]), bb.max()[d]);
            }
        }
    }

    if (writeFields_)
    {
        forAll(phaseNames_, phasei)
        {
            const word fieldName
            (
                IOobject::groupName("rCFDdisplacement", phaseNames_[phasei])
            );

            if (!mesh_.foundObject<volVectorField>(fieldName))
            {
                volVectorField* dPtr = new volVectorField
                (
                    IOobject
                    (
                        fieldName,
                        time_.name(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector(dimLength, Zero),
                    zeroGradientFvPatchVectorField::typeName
                );

                dPtr->store();
            }

            volVectorField& displacement =
                mesh_.lookupObjectRef<volVectorField>(fieldName);

            for (direction d = 0; d < 3; d++)
            {
                displacement.primitiveFieldRef().replace
                (
                    d,
                    mesh_.C().primitiveField().component(d)
                  - labels_[3*phasei + d].primitiveField()
                );
            }

            displacement.correctBoundaryConditions();
        }
    }

    rcfd_of_prep_frame(vof.begin(), origin.begin());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::rCFDPrep::rCFDPrep
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    recordStart_(0),
    userLibrary_("rCFD/librCFD_user_prep.so"),
    normPhase_(-1),
    frameInterval_(0),
    nFrames_(0),
    residualAlpha_(1e-6),
    labelMaxCo_(1),
    maxLabelSubCycles_(200),
    coFrame_(0),
    nSubFrame_(1),
    writeFields_(false),
    labels_(),
    opened_(false),
    recording_(false),
    done_(false),
    nextFrameTime_(great),
    nLateFrames_(0),
    cpuLabels_(0),
    coronaFace_()
{
    if (functionObject::postProcess)
    {
        FatalErrorInFunction
            << type() << " works during a solver run only"
            << exit(FatalError);
    }

    read(dict);

    const word cmptNames[3] = {"x", "y", "z"};

    labels_.setSize(3*phaseNames_.size());

    forAll(phaseNames_, phasei)
    {
        for (direction d = 0; d < 3; d++)
        {
            labels_.set
            (
                3*phasei + d,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            word("rCFD" + cmptNames[d]),
                            phaseNames_[phasei]
                        ),
                        time_.name(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar(dimLength, 0),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
        }
    }

    resetLabels();

    // Mesh to the library; the user header of the case sets the frames
    open();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::rCFDPrep::~rCFDPrep()
{
    // OpenFOAM does not call end() when a run is stopped by a change of
    // controlDict (endTime, stopAt) while it runs: close the database here,
    // so that the frames written so far can be used
    if (opened_ && !done_)
    {
        rcfd_of_prep_close();

        done_ = true;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::rCFDPrep::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if
    (
        dict.found("executeControl")
     || dict.found("writeControl")
     || dict.found("outputControl")
    )
    {
        FatalIOErrorInFunction(dict)
            << "Remove executeControl/writeControl from the dictionary of "
            << name() << ": the labels have to be advanced every time step"
            << exit(FatalIOError);
    }

    if (opened_)
    {
        // The library is open: the settings of the database cannot change
        return true;
    }

    if (dict.found("frameInterval") || dict.found("nFrames"))
    {
        FatalIOErrorInFunction(dict)
            << "frameInterval and nFrames are set in the user header "
            << "rCFD/rCFD_user_prep.h (global_time_step, number_of_frames), "
            << "not here. Remove them from the dictionary of " << name()
            << exit(FatalIOError);
    }

    if (dict.found("phases"))
    {
        phaseNames_ = dict.lookup<wordList>("phases");
    }
    else
    {
        phaseNames_.setSize(2);
        phaseNames_[0] = dict.lookup<word>("gasPhase");
        phaseNames_[1] = dict.lookup<word>("solidPhase");
    }

    if (phaseNames_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "No phases given" << exit(FatalIOError);
    }

    const dictionary& fieldsDict = dict.subOrEmptyDict("fields");

    alphaNames_.setSize(phaseNames_.size());
    fluxNames_.setSize(phaseNames_.size());
    rhoNames_.setSize(phaseNames_.size());

    forAll(phaseNames_, phasei)
    {
        const word& phase = phaseNames_[phasei];
        const dictionary& phaseDict = fieldsDict.subOrEmptyDict(phase);

        alphaNames_[phasei] =
            phaseDict.lookupOrDefault<word>
            (
                "alpha",
                IOobject::groupName("alpha", phase)
            );

        fluxNames_[phasei] =
            phaseDict.lookupOrDefault<word>
            (
                "flux",
                IOobject::groupName("alphaRhoPhi", phase)
            );

        rhoNames_[phasei] =
            phaseDict.lookupOrDefault<word>
            (
                "rho",
                IOobject::groupName("rho", phase)
            );
    }

    recordStart_ = dict.lookup<scalar>("recordStart");

    userLibrary_ =
        dict.lookupOrDefault<fileName>
        (
            "userLibrary",
            "rCFD/librCFD_user_prep.so"
        );

    normPhase_ = dict.lookupOrDefault<label>("normPhase", -1);

    divScheme_ =
        dict.lookupOrDefault<string>("divScheme", "Gauss limitedLinear 1");

    solverDict_ = dict.subOrEmptyDict("solver");

    if (!solverDict_.found("solver"))
    {
        solverDict_.add("solver", word("PBiCGStab"));
        solverDict_.add("preconditioner", word("DILU"));
    }
    if (!solverDict_.found("tolerance"))
    {
        solverDict_.add("tolerance", scalar(1e-10));
    }
    if (!solverDict_.found("relTol"))
    {
        solverDict_.add("relTol", scalar(0));
    }

    residualAlpha_ = dict.lookupOrDefault<scalar>("residualAlpha", 1e-6);
    labelMaxCo_ = dict.lookupOrDefault<scalar>("labelMaxCo", 1);
    maxLabelSubCycles_ = dict.lookupOrDefault<label>("maxLabelSubCycles", 200);

    if (labelMaxCo_ <= 0 || maxLabelSubCycles_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "labelMaxCo must be positive and maxLabelSubCycles at least 1"
            << exit(FatalIOError);
    }
    writeFields_ = dict.lookupOrDefault<Switch>("writeFields", false);

    return true;
}


Foam::wordList Foam::functionObjects::rCFDPrep::fields() const
{
    DynamicList<word> names;

    forAll(phaseNames_, phasei)
    {
        if (alphaNames_[phasei] != "none")
        {
            names.append(alphaNames_[phasei]);
        }

        names.append(fluxNames_[phasei]);
    }

    return wordList(names);
}


bool Foam::functionObjects::rCFDPrep::execute()
{
    if (done_)
    {
        return true;
    }

    const scalar t = time_.value();
    const scalar tol = 1e-6*frameInterval_;

    if (!recording_)
    {
        const scalar tOld = t - time_.deltaTValue();

        if
        (
            time_.timeIndex() == time_.startTimeIndex() + 1
         && tOld >= recordStart_ - tol
        )
        {
            // The run started at or after recordStart (restart after the
            // spin-up): the labels are at the cell centres since the start,
            // so the first time step belongs to the first frame
            recording_ = true;
            nextFrameTime_ = tOld + frameInterval_;

            Info<< type() << " " << name() << ": recording starts at t = "
                << tOld << " s (start of the run)" << endl;
        }
        else
        {
            if (t >= recordStart_ - tol)
            {
                resetLabels();

                recording_ = true;
                nextFrameTime_ = t + frameInterval_;

                Info<< type() << " " << name() << ": recording starts at t = "
                    << t << " s" << endl;
            }

            return true;
        }
    }

    // Advance the labels over the time step just solved
    {
        cpuTime timer;

        const int oldDebug = solverPerformance::debug;
        solverPerformance::debug = 0;

        forAll(phaseNames_, phasei)
        {
            solveLabels(phasei);
        }

        solverPerformance::debug = oldDebug;

        cpuLabels_ += timer.elapsedCpuTime();
    }

    if (t >= nextFrameTime_ - tol)
    {
        if (t > nextFrameTime_ + tol)
        {
            nLateFrames_++;

            if (nLateFrames_ <= 5)
            {
                WarningInFunction
                    << "Frame taken at t = " << t << " s, "
                    << t - nextFrameTime_ << " s after its time. Use "
                    << "adjustTimeStep yes or a time step that divides the "
                    << "frame interval." << endl;
            }
        }

        Info<< type() << " " << name() << ": frame "
            << rcfd_of_prep_frames_written() << ": label Courant number up to "
            << coFrame_ << ", up to " << nSubFrame_
            << " label sub-steps per time step"
            << (
                   coFrame_ > labelMaxCo_*maxLabelSubCycles_*(1 + 1e-6)
                 ? " (limited by maxLabelSubCycles)"
                 : ""
               )
            << endl;

        coFrame_ = 0;
        nSubFrame_ = 1;

        handOverFrame();
        resetLabels();

        nextFrameTime_ += frameInterval_;

        if (rcfd_of_prep_frames_written() >= nFrames_)
        {
            rcfd_of_prep_close();

            done_ = true;

            Info<< type() << " " << name() << ": database complete at t = "
                << t << " s (" << nFrames_ << " frames, label transport "
                << cpuLabels_ << " s cpu on the master). The simulation can "
                << "be stopped." << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::rCFDPrep::write()
{
    // Displacement of the last complete frame, written with the solver fields
    if (writeFields_ && time_.writeTime())
    {
        forAll(phaseNames_, phasei)
        {
            const word fieldName
            (
                IOobject::groupName("rCFDdisplacement", phaseNames_[phasei])
            );

            if (mesh_.foundObject<volVectorField>(fieldName))
            {
                mesh_.lookupObject<volVectorField>(fieldName).write();
            }
        }
    }

    return true;
}


bool Foam::functionObjects::rCFDPrep::end()
{
    if (opened_ && !done_)
    {
        rcfd_of_prep_close();

        done_ = true;
    }

    return true;
}


Foam::scalar Foam::functionObjects::rCFDPrep::timeToNextAction()
{
    if (done_)
    {
        return vGreat;
    }

    const scalar t = time_.value();

    if (!recording_)
    {
        return (recordStart_ > t) ? recordStart_ - t : vGreat;
    }

    return max(nextFrameTime_ - t, scalar(0));
}


// ************************************************************************* //
