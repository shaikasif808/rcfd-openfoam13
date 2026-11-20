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
#include "cyclicPolyPatch.H"
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "cyclicFvPatch.H"
#include "HashTable.H"
#include "labelPair.H"
#include "syncTools.H"
#include "labelIOList.H"
#include "IStringStream.H"
#include "cpuTime.H"
#include "OSspecific.H"

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
    neighbourValues(cellValues, nbrValues);

    const label nImg = imageFace_.size();

    forAll(cellValues, celli)
    {
        values[stride*celli + offset] = cellValues[celli];
    }

    // An image cell carries the value of the cell on the other side
    forAll(imageFace_, i)
    {
        values[stride*(nInt + i) + offset] = nbrValues[imageFace_[i]];
    }

    forAll(coronaFace_, i)
    {
        values[stride*(nInt + nImg + i) + offset] = nbrValues[coronaFace_[i]];
    }
}


void Foam::functionObjects::rCFDPrep::setBoundaryFaces()
{
    const fvBoundaryMesh& fvb = mesh_.boundary();

    bfStart_.setSize(fvb.size() + 1);

    label n = 0;

    forAll(fvb, patchi)
    {
        bfStart_[patchi] = n;
        n += fvb[patchi].size();
    }

    bfStart_[fvb.size()] = n;

    bfCell_.setSize(n);
    bfSf_.setSize(n);

    forAll(fvb, patchi)
    {
        const labelUList& faceCells = fvb[patchi].faceCells();
        const vectorField& SfP = mesh_.Sf().boundaryField()[patchi];

        forAll(faceCells, i)
        {
            bfCell_[bfStart_[patchi] + i] = faceCells[i];
            bfSf_[bfStart_[patchi] + i] = SfP[i];
        }
    }
}


void Foam::functionObjects::rCFDPrep::neighbourValues
(
    const scalarField& cellValues,
    scalarList& nbrValues
) const
{
    // A field with the coupled conditions of the mesh: its coupled patches
    // know the value on the other side, whatever the kind of coupling
    volScalarField f
    (
        IOobject
        (
            "rCFDneighbour",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    );

    f.primitiveFieldRef() = cellValues;
    f.correctBoundaryConditions();

    nbrValues.setSize(bfCell_.size());
    nbrValues = 0;

    forAll(f.boundaryField(), patchi)
    {
        if (!mesh_.boundary()[patchi].coupled()) continue;

        const scalarField nbr(f.boundaryField()[patchi].patchNeighbourField());

        forAll(nbr, i)
        {
            nbrValues[bfStart_[patchi] + i] = nbr[i];
        }
    }
}


void Foam::functionObjects::rCFDPrep::setPeriodic()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Area vector and area-weighted centre of every cyclic patch, over all
    // processors. The faces of a cyclic patch that lie on a processor
    // boundary are in processorCyclic patches.
    // The patches that are not processor patches come first and are the
    // same on all processors: lists of that length can be combined
    label nGlobalPatches = 0;

    forAll(patches, patchi)
    {
        if (!isA<processorPolyPatch>(patches[patchi]))
        {
            nGlobalPatches = patchi + 1;
        }
    }

    vectorField sumSf(nGlobalPatches, Zero), sumC(nGlobalPatches, Zero);
    scalarField sumA(nGlobalPatches, 0);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        label cyci = -1;

        if (isA<processorCyclicPolyPatch>(pp))
        {
            cyci = refCast<const processorCyclicPolyPatch>(pp).referPatchIndex();
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            cyci = patchi;
        }

        if (cyci < 0) continue;

        // Faces of the finite volume patch: for a non-conformal patch the
        // coupled faces
        const vectorField& Sf = mesh_.Sf().boundaryField()[patchi];
        const vectorField& Cf = mesh_.Cf().boundaryField()[patchi];

        forAll(Sf, i)
        {
            sumSf[cyci] += Sf[i];
            sumC[cyci] += mag(Sf[i])*Cf[i];
            sumA[cyci] += mag(Sf[i]);
        }
    }

    Pstream::listCombineGather(sumSf, plusEqOp<vector>());
    Pstream::listCombineScatter(sumSf);
    Pstream::listCombineGather(sumC, plusEqOp<vector>());
    Pstream::listCombineScatter(sumC);
    Pstream::listCombineGather(sumA, plusEqOp<scalar>());
    Pstream::listCombineScatter(sumA);

    DynamicList<point> centres;
    DynamicList<vector> normals, shifts;
    string names;

    periodicPlaneOfPatch_.setSize(patches.size());
    periodicPlaneOfPatch_ = -1;

    interfacePatch_.setSize(patches.size());
    interfacePatch_ = false;

    string interfaceNames;

    forAll(patches, patchi)
    {
        if (!isA<cyclicPolyPatch>(patches[patchi])) continue;

        const cyclicPolyPatch& cpp =
            refCast<const cyclicPolyPatch>(patches[patchi]);

        const label nbri = cpp.nbrPatchIndex();

        if (sumA[patchi] <= 0 || sumA[nbri] <= 0) continue;

        if (cpp.transform().transforms())
        {
            FatalErrorInFunction
                << "Patch " << cpp.name() << ": periodic patches with a "
                << "rotation are not supported" << exit(FatalError);
        }

        const vector shift
        (
            sumC[nbri]/sumA[nbri] - sumC[patchi]/sumA[patchi]
        );

        // No shift: an interface inside the domain, not a periodic side
        if (mag(shift) < 1e-6*sqrt(sumA[patchi]))
        {
            interfacePatch_[patchi] = true;

            interfaceNames +=
                (interfaceNames.size() ? " " : "") + cpp.name();

            continue;
        }

        // Flat and opposite to its neighbour
        if
        (
            mag(sumSf[patchi]) < 0.9999*sumA[patchi]
         || mag(sumSf[patchi] + sumSf[nbri]) > 1e-4*sumA[patchi]
        )
        {
            FatalErrorInFunction
                << "Patch " << cpp.name() << ": a periodic patch must be "
                << "plane and parallel to its neighbour patch "
                << patches[nbri].name() << exit(FatalError);
        }

        periodicPlaneOfPatch_[patchi] = centres.size();

        centres.append(sumC[patchi]/sumA[patchi]);
        normals.append(sumSf[patchi]/mag(sumSf[patchi]));
        shifts.append(shift);

        names += (names.size() ? " " : "") + cpp.name();
    }

    periodicCentre_ = centres;
    periodicNormal_ = normals;
    periodicShift_ = shifts;

    displacementForm_ =
        labelForm_ == "displacement"
     || (labelForm_ == "auto" && periodicShift_.size());

    if (periodicShift_.size() && !displacementForm_)
    {
        FatalErrorInFunction
            << "The mesh has periodic patches (" << names.c_str() << "): "
            << "labelForm position cannot be used, the position jumps "
            << "across them" << exit(FatalError);
    }

    if (interfaceNames.size())
    {
        Info<< type() << " " << name() << ": interfaces inside the domain "
            << interfaceNames.c_str() << ": the cells on the two sides are "
            << "neighbours" << nl;
    }

    if (periodicShift_.size())
    {
        Info<< type() << " " << name() << ": periodic patches "
            << names.c_str() << nl;

        forAll(periodicShift_, k)
        {
            Info<< "    plane through " << periodicCentre_[k] << ", normal "
                << periodicNormal_[k] << ": origin points beyond it are "
                << "moved by " << periodicShift_[k] << nl;
        }
    }

    Info<< type() << " " << name() << ": labels in "
        << (displacementForm_ ? "displacement" : "position") << " form"
        << endl;
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

    // Values of the cells on the other side of the faces of the coupled
    // patches (index: list of all patch faces, see setBoundaryFaces)
    labelList nbrGid;
    {
        scalarField gidS(nInt);
        forAll(gidS, celli)
        {
            gidS[celli] = gid[celli];
        }

        scalarList nbr;
        neighbourValues(gidS, nbr);

        nbrGid.setSize(nbr.size());
        forAll(nbr, i)
        {
            nbrGid[i] = label(nbr[i] + 0.5);
        }
    }

    const vectorField& C = mesh_.cellCentres();
    const scalarField& V = mesh_.cellVolumes();
    const vectorField& Sf = mesh_.faceAreas();

    scalarList nbrV;
    neighbourValues(V, nbrV);

    FixedList<scalarList, 3> nbrC;
    for (direction d = 0; d < 3; d++)
    {
        const scalarField Cd(C.component(d));
        neighbourValues(Cd, nbrC[d]);
    }

    const fvBoundaryMesh& fvb = mesh_.boundary();

    // Partition faces: the faces of the processor patches
    DynamicList<label> pfFace;
    DynamicList<label> pfHost;

    // Faces of interfaces inside the domain with both cells on this
    // processor
    DynamicList<label> ifFace;

    forAll(fvb, patchi)
    {
        const fvPatch& fvp = fvb[patchi];

        if (isA<processorCyclicFvPatch>(fvp) || isA<cyclicFvPatch>(fvp))
        {
            const bool onProcessor = isA<processorCyclicFvPatch>(fvp);

            const label cyci =
                onProcessor
              ? refCast<const processorCyclicPolyPatch>
                (
                    fvp.patch()
                ).referPatchIndex()
              : patchi;

            if (!interfacePatch_[cyci])
            {
                // Periodic (conformal or not): no face connection for the
                // library here. The transport across these patches is in
                // the labels and the image cells, see setPeriodic()
            }
            else if (onProcessor)
            {
                const label host =
                    refCast<const processorFvPatch>(fvp).neighbProcNo();

                forAll(fvp, i)
                {
                    pfFace.append(bfStart_[patchi] + i);
                    pfHost.append(host);
                }
            }
            else
            {
                forAll(fvp, i)
                {
                    ifFace.append(bfStart_[patchi] + i);
                }
            }
        }
        else if (isA<processorFvPatch>(fvp))
        {
            const label host =
                refCast<const processorFvPatch>(fvp).neighbProcNo();

            forAll(fvp, i)
            {
                pfFace.append(bfStart_[patchi] + i);
                pfHost.append(host);
            }
        }
        else if (fvp.coupled())
        {
            FatalErrorInFunction
                << "Patch " << fvp.name() << " of type " << fvp.type()
                << ": coupled patches other than processor and cyclic "
                << "patches are not supported" << exit(FatalError);
        }
    }

    // Corona cells: the distinct cells behind the partition faces, ordered
    // by host and then by global cell index
    std::vector<std::pair<long long, label>> keys(pfFace.size());

    forAll(pfFace, i)
    {
        keys[i].first =
            static_cast<long long>(pfHost[i])*static_cast<long long>(nGlobal)
          + static_cast<long long>(nbrGid[pfFace[i]]);
        keys[i].second = i;
    }

    std::sort(keys.begin(), keys.end());

    labelList pfCorona(pfFace.size(), -1);
    DynamicList<label> coronaFace;

    for (std::size_t k = 0; k < keys.size(); k++)
    {
        if (k == 0 || keys[k].first != keys[k - 1].first)
        {
            coronaFace.append(pfFace[keys[k].second]);
        }

        pfCorona[keys[k].second] = coronaFace.size() - 1;
    }

    coronaFace_ = coronaFace;

    // Faces across the interfaces inside the domain: one per pair of
    // cells, from the side of the cell with the lower index
    DynamicList<label> ifCell0, ifCell1;
    DynamicList<vector> ifSf;

    if (ifFace.size())
    {
        Map<label> localCell(2*nInt);

        forAll(gid, celli)
        {
            localCell.insert(gid[celli], celli);
        }

        HashTable<label, labelPair, labelPair::Hash<>> faceOf;

        forAll(ifFace, i)
        {
            const label bFacei = ifFace[i];
            const labelPair cells(bfCell_[bFacei], localCell[nbrGid[bFacei]]);

            if (cells.first() >= cells.second()) continue;

            if (faceOf.found(cells))
            {
                ifSf[faceOf[cells]] += bfSf_[bFacei];
            }
            else
            {
                faceOf.insert(cells, ifCell0.size());

                ifCell0.append(cells.first());
                ifCell1.append(cells.second());
                ifSf.append(bfSf_[bFacei]);
            }
        }
    }

    // Image cells: behind every periodic patch one copy of every cell on
    // the other side that shares a face with a cell here. Across a
    // non-conformal patch a cell can face several cells: one image per
    // cell on the other side, one face per pair of cell and image.
    DynamicList<label> imgFace, imgHost, imgPlane;
    DynamicList<label> imgFaceCell, imgFaceImage;
    DynamicList<vector> imgFaceSf;

    if (periodicImages_ && periodicShift_.size())
    {
        // (plane, cell on the other side) -> image, (cell, image) -> face
        HashTable<label, labelPair, labelPair::Hash<>> imageOf, faceOf;

        forAll(fvb, patchi)
        {
            const fvPatch& fvp = fvb[patchi];
            const polyPatch& pp = fvp.patch();

            label cyci = -1;
            label hostProc = Pstream::myProcNo();

            if (isA<processorCyclicPolyPatch>(pp))
            {
                const processorCyclicPolyPatch& pcpp =
                    refCast<const processorCyclicPolyPatch>(pp);

                cyci = pcpp.referPatchIndex();
                hostProc = pcpp.neighbProcNo();
            }
            else if (isA<cyclicPolyPatch>(pp))
            {
                cyci = patchi;
            }

            if (cyci < 0 || periodicPlaneOfPatch_[cyci] < 0) continue;

            const label planei = periodicPlaneOfPatch_[cyci];

            forAll(fvp, i)
            {
                const label bFacei = bfStart_[patchi] + i;

                const labelPair imageKey(planei, nbrGid[bFacei]);

                label imagei = -1;

                if (imageOf.found(imageKey))
                {
                    imagei = imageOf[imageKey];
                }
                else
                {
                    imagei = imgFace.size();
                    imageOf.insert(imageKey, imagei);

                    imgFace.append(bFacei);
                    imgHost.append(hostProc);
                    imgPlane.append(planei);
                }

                const labelPair faceKey(bfCell_[bFacei], imagei);

                if (faceOf.found(faceKey))
                {
                    imgFaceSf[faceOf[faceKey]] += bfSf_[bFacei];
                }
                else
                {
                    faceOf.insert(faceKey, imgFaceCell.size());

                    imgFaceCell.append(bfCell_[bFacei]);
                    imgFaceImage.append(imagei);
                    imgFaceSf.append(bfSf_[bFacei]);
                }
            }
        }
    }

    imageFace_ = imgFace;

    const label nImg = imageFace_.size();
    const label nImgFaces = imgFaceCell.size();
    const label nExt = coronaFace_.size();
    const label nCells = nInt + nImg + nExt;
    // Faces between cells of this processor: mesh faces and interfaces
    const label nOwnFaces = nIntFaces + ifCell0.size();
    const label nFaces = nOwnFaces + nImgFaces + pfFace.size();

    // Arrays for the library: cells of this processor, image cells, corona
    // cells; faces between own cells, faces to image cells, partition faces
    std::vector<double> x(3*nCells), vol(nCells), area(3*nFaces);
    std::vector<int> host(nCells), gidAll(nCells), c0(nFaces), c1(nFaces);
    std::vector<int> imageHost(nImg > 0 ? nImg : 1);

    imageCentre_.setSize(nImg);
    periodicDepth_.setSize(periodicShift_.size());
    periodicDepth_ = 0;

    forAll(imageFace_, i)
    {
        const label bFacei = imageFace_[i];
        const label k = imgPlane[i];
        const label cellj = nInt + i;

        // Where the cell on the other side would be if the mesh continued
        const point p
        (
            point(nbrC[0][bFacei], nbrC[1][bFacei], nbrC[2][bFacei])
          - periodicShift_[k]
        );

        imageCentre_[i] = p;

        for (direction d = 0; d < 3; d++)
        {
            x[3*cellj + d] = p[d];
        }

        vol[cellj] = nbrV[bFacei];
        host[cellj] = Pstream::myProcNo();
        gidAll[cellj] = nbrGid[bFacei];
        imageHost[i] = imgHost[i];

        periodicDepth_[k] = max
        (
            periodicDepth_[k],
            2*((p - periodicCentre_[k]) & periodicNormal_[k])
        );
    }

    Pstream::listCombineGather(periodicDepth_, maxEqOp<scalar>());
    Pstream::listCombineScatter(periodicDepth_);

    if (periodicShift_.size())
    {
        Info<< type() << " " << name() << ": "
            << returnReduce(nImg, sumOp<label>())
            << " image cells behind the periodic patches, depth of the "
            << "image layers " << periodicDepth_ << endl;
    }

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
        const label bFacei = pfFace[i];
        const label cellj = nInt + nImg + pfCorona[i];

        // Every face of a corona cell carries the same values
        for (direction d = 0; d < 3; d++)
        {
            x[3*cellj + d] = nbrC[d][bFacei];
        }

        vol[cellj] = nbrV[bFacei];
        host[cellj] = pfHost[i];
        gidAll[cellj] = nbrGid[bFacei];
    }

    forAll(ifCell0, i)
    {
        const label k = nIntFaces + i;

        c0[k] = ifCell0[i];
        c1[k] = ifCell1[i];

        for (direction d = 0; d < 3; d++)
        {
            area[3*k + d] = ifSf[i][d];
        }
    }

    forAll(imgFaceCell, i)
    {
        const label k = nOwnFaces + i;

        c0[k] = imgFaceCell[i];
        c1[k] = nInt + imgFaceImage[i];

        for (direction d = 0; d < 3; d++)
        {
            area[3*k + d] = imgFaceSf[i][d];
        }
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
        const label bFacei = pfFace[i];
        const label k = nOwnFaces + nImgFaces + i;

        c0[k] = bfCell_[bFacei];
        c1[k] = nInt + nImg + pfCorona[i];

        // Area vector of a boundary face points out of the own cell
        for (direction d = 0; d < 3; d++)
        {
            area[3*k + d] = bfSf_[bFacei][d];
        }
    }

    rcfd_of_mesh_t m;

    m.n_int = nInt + nImg;
    m.n_ext = nExt;
    m.n_fint = nOwnFaces + nImgFaces;
    m.n_fext = pfFace.size();
    m.n_img = nImg;
    m.n_fimg = nImgFaces;
    m.img_host = imageHost.data();
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

    List<double> extra(extraIntervals_.size());
    forAll(extra, i)
    {
        extra[i] = extraIntervals_[i];
    }

    s.n_extra = extra.size();
    s.extra_intervals = extra.size() ? extra.begin() : nullptr;

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

        // Faces of the finite volume patches: where a patch is covered by
        // a non-conformal couple its face areas are (nearly) zero
        const fvBoundaryMesh& fvb = mesh_.boundary();

        label nPatches = 0;
        forAll(fvb, patchi)
        {
            if (!fvb[patchi].coupled()) nPatches++;
        }

        os << nPatches << '\n';

        forAll(fvb, patchi)
        {
            const fvPatch& fvp = fvb[patchi];

            if (fvp.coupled()) continue;

            // The error patches of non-conformal couples have no faces in
            // the polyhedral mesh: they are listed without faces
            const bool faces = fvp.size() == fvp.patch().size();

            const labelUList& faceCells = fvp.faceCells();
            const scalarField& magSf = mesh_.magSf().boundaryField()[patchi];
            const vectorField& Cf = mesh_.Cf().boundaryField()[patchi];

            if (!faces)
            {
                // An empty patch has faces in the polyhedral mesh only
                const polyPatch& pp = fvp.patch();
                const labelUList& ppCells = pp.faceCells();
                const vectorField ppSf(pp.faceAreas());
                const vectorField ppCf(pp.faceCentres());

                os << pp.name() << ' ' << pp.type() << ' ' << pp.size()
                   << '\n';

                forAll(pp, i)
                {
                    os << ppCells[i] << ' ' << mag(ppSf[i]) << ' '
                       << ppCf[i].x() << ' ' << ppCf[i].y() << ' '
                       << ppCf[i].z() << '\n';
                }

                continue;
            }

            os << fvp.name() << ' ' << fvp.patch().type() << ' '
               << fvp.size() << '\n';

            forAll(fvp, i)
            {
                os << faceCells[i] << ' ' << magSf[i] << ' '
                   << Cf[i].x() << ' ' << Cf[i].y() << ' '
                   << Cf[i].z() << '\n';
            }
        }
    }

    frameInterval_ = info.frame_interval;
    nFrames_ = info.n_frames;

    Info<< type() << " " << name() << ": " << nFrames_ << " frames of "
        << frameInterval_ << " s from t = " << recordStart_ << " s"
        << " (rCFD_user_prep.h)" << endl;

    forAll(setInterval_, seti)
    {
        setInterval_[seti] = info.db_interval[seti];
        setFrames_[seti] = info.db_frames[seti];
        setDir_[seti] = fileName(info.db_dir[seti]);

        if (seti > 0)
        {
            // The patch file of the run mode belongs to every database
            const word patchName
            (
                "patches_" + Foam::name(Pstream::myProcNo())
            );

            Foam::cp
            (
                fileName("data/topo")/patchName,
                setDir_[seti]/"data/topo"/patchName
            );

            Info<< type() << " " << name() << ": " << setFrames_[seti]
                << " frames of " << setInterval_[seti] << " s in "
                << setDir_[seti].c_str() << "data (extraFrameIntervals)"
                << endl;
        }
    }

    opened_ = true;
}


void Foam::functionObjects::rCFDPrep::resetLabels(const label seti)
{
    const label nPerSet = 3*phaseNames_.size();

    forAll(labels_, i)
    {
        if (seti >= 0 && i/nPerSet != seti) continue;

        volScalarField& X = labels_[i];

        if (displacementForm_)
        {
            X.primitiveFieldRef() = 0;
        }
        else
        {
            X.primitiveFieldRef() =
                mesh_.C().primitiveField().component(i % 3);
        }

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
        neighbourValues(m, nbrM);

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];
            const labelUList& faceCells = fvp.faceCells();
            const scalarField& fluxP = flux.boundaryField()[patchi];

            forAll(fluxP, i)
            {
                const label celli = faceCells[i];

                const scalar mf =
                    fvp.coupled()
                  ? max(m[celli], nbrM[bfStart_[patchi] + i])
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

        forAll(coFrame_, seti)
        {
            coFrame_[seti] = max(coFrame_[seti], co);
            nSubFrame_[seti] = max(nSubFrame_[seti], nSub);
        }
    }

    // The displacement form always goes through the sub-step form of the
    // equation (one sub-step is the plain time step)
    if (nSub > 1 || displacementForm_)
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

        // Displacement form: with X = x + xi the label equation becomes
        //   ddt(m xi) + div(flux xi) - contErr xi
        //       = - 1/V sum_faces flux_f (x_f - x_cell)
        // The right-hand side holds only vectors inside a cell, so xi is
        // continuous across periodic patches. Faces of boundary patches
        // that are not coupled give nothing: this is the zero-gradient
        // condition of the position form (the material that enters
        // through an inlet comes from the cell itself).
        tmp<DimensionedField<vector, volMesh>> tSxi;

        if (displacementForm_)
        {
            tSxi = new DimensionedField<vector, volMesh>
            (
                IOobject("Sxi", time_.name(), mesh_),
                mesh_,
                dimensionedVector(flux.dimensions()*dimLength/dimVolume, Zero)
            );

            vectorField& S = tSxi.ref().primitiveFieldRef();

            const vectorField& C = mesh_.C().primitiveField();
            const vectorField& Cf = mesh_.Cf().primitiveField();
            const scalarField& fluxI = flux.primitiveField();
            const labelUList& own = mesh_.owner();
            const labelUList& nei = mesh_.neighbour();

            forAll(nei, facei)
            {
                S[own[facei]] -= fluxI[facei]*(Cf[facei] - C[own[facei]]);
                S[nei[facei]] += fluxI[facei]*(Cf[facei] - C[nei[facei]]);
            }

            forAll(mesh_.boundary(), patchi)
            {
                const fvPatch& fvp = mesh_.boundary()[patchi];

                if (!fvp.coupled()) continue;

                const labelUList& faceCells = fvp.faceCells();
                const scalarField& fluxP = flux.boundaryField()[patchi];
                const vectorField& CfP = mesh_.Cf().boundaryField()[patchi];

                forAll(faceCells, i)
                {
                    S[faceCells[i]] -= fluxP[i]*(CfP[i] - C[faceCells[i]]);
                }
            }

            S /= mesh_.V();
        }

        for (label subi = 0; subi < nSub; subi++)
        {
            const scalar f0 = scalar(subi)/nSub;
            const scalar f1 = scalar(subi + 1)/nSub;

            const dField c0("c0", (mOld + f0*(mNew - mOld) + mRes)*rDtSub);
            const dField c1("c1", (mOld + f1*(mNew - mOld) + mRes)*rDtSub);

            // Every label set that is still recording, with the same flux
            forAll(setDone_, seti)
            {
                if (setDone_[seti]) continue;

                for (direction d = 0; d < 3; d++)
                {
                    if (solutionD[d] == -1) continue;

                    volScalarField& X = labels_[labelIndex(seti, phasei, d)];

                    IStringStream scheme(divScheme_);

                    tmp<fvScalarMatrix> tXEqn
                    (
                        fv::convectionScheme<scalar>::New
                        (
                            mesh_,
                            flux,
                            scheme
                        )().fvmDiv(flux, X)
                      + fvm::Sp(c1 - contErr.internalField(), X)
                      - fvm::Su(c0*X.internalField(), X)
                    );

                    if (displacementForm_)
                    {
                        tXEqn.ref() -= fvm::Su(tSxi().component(d), X);
                    }

                    tXEqn.ref().solve(solverDict_);
                }
            }
        }

        return;
    }

    forAll(setDone_, seti)
    {
      if (setDone_[seti]) continue;

      for (direction d = 0; d < 3; d++)
      {
        if (solutionD[d] == -1) continue;

        volScalarField& X = labels_[labelIndex(seti, phasei, d)];

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
}


void Foam::functionObjects::rCFDPrep::handOverFrame(const label seti)
{
    const label nInt = mesh_.nCells();
    const label nImg = imageFace_.size();
    const label nCells = nInt + nImg + coronaFace_.size();
    const label nPhases = phaseNames_.size();

    // vof[cell*nPhases + phase], origin[(phase*nCells + cell)*3 + direction]
    List<double> vof(nCells*nPhases), origin(nPhases*3*nCells);
    List<double> displ(nPhases*3*nCells);

    // Box of the whole mesh
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
            if (displacementForm_)
            {
                // Position the material came from; may lie beyond a
                // periodic patch
                const scalarField Xd
                (
                    mesh_.C().primitiveField().component(d)
                  + labels_[labelIndex(seti, phasei, d)].primitiveField()
                );

                cellAndCoronaValues(Xd, origin, 3, offset + d);
            }
            else
            {
                cellAndCoronaValues
                (
                    labels_[labelIndex(seti, phasei, d)].primitiveField(),
                    origin,
                    3,
                    offset + d
                );
            }
        }

        // Displacement before the origin is moved through a periodic patch
        // (for the report of the bridge)
        for (label i = 0; i < nCells; i++)
        {
            for (direction d = 0; d < 3; d++)
            {
                displ[offset + 3*i + d] = origin[offset + 3*i + d];
            }
        }

        // An image cell is its own origin
        forAll(imageFace_, i)
        {
            for (direction d = 0; d < 3; d++)
            {
                origin[offset + 3*(nInt + i) + d] = imageCentre_[i][d];
                displ[offset + 3*(nInt + i) + d] = imageCentre_[i][d];
            }
        }

        // An origin beyond a periodic patch is on the other side of the
        // domain. There the library finds a cell with all its neighbours;
        // an image cell has only the cells in front of it as neighbours
        // and is for the interpolation, not a source to start from.
        List<bool> inImageLayer(nCells, false);

        if (periodicShift_.size())
        {
            for (label i = 0; i < nCells; i++)
            {
                if (i >= nInt && i < nInt + nImg) continue;

                point p
                (
                    origin[offset + 3*i],
                    origin[offset + 3*i + 1],
                    origin[offset + 3*i + 2]
                );

                for (label pass = 0; pass < 2; pass++)
                {
                    forAll(periodicShift_, k)
                    {
                        const scalar beyond =
                            (p - periodicCentre_[k]) & periodicNormal_[k];

                        if (beyond > 0)
                        {
                            p += periodicShift_[k];
                        }
                    }
                }

                forAll(periodicShift_, k)
                {
                    if (((p - periodicCentre_[k]) & periodicNormal_[k]) > 0)
                    {
                        inImageLayer[i] = true;
                    }
                }

                origin[offset + 3*i] = p.x();
                origin[offset + 3*i + 1] = p.y();
                origin[offset + 3*i + 2] = p.z();
            }
        }

        // An origin cannot lie outside the domain (a limited scheme may
        // overshoot by a small amount next to a wall)
        for (label i = 0; i < nCells; i++)
        {
            if (inImageLayer[i] || (i >= nInt && i < nInt + nImg)) continue;

            for (direction d = 0; d < 3; d++)
            {
                double& x = origin[offset + 3*i + d];
                x = min(max(x, bb.min()[d]), bb.max()[d]);
            }
        }
    }

    // Displacement fields for ParaView: those of the database of the header
    if (writeFields_ && seti == 0)
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
                if (displacementForm_)
                {
                    displacement.primitiveFieldRef().replace
                    (
                        d,
                        -labels_[3*phasei + d].primitiveField()
                    );
                }
                else
                {
                    displacement.primitiveFieldRef().replace
                    (
                        d,
                        mesh_.C().primitiveField().component(d)
                      - labels_[3*phasei + d].primitiveField()
                    );
                }
            }

            displacement.correctBoundaryConditions();
        }
    }

    // displ: unwrapped origin minus the position of the cell or corona cell
    {
        scalarList nbrC;

        for (direction d = 0; d < 3; d++)
        {
            const scalarField Cd(mesh_.C().primitiveField().component(d));
            neighbourValues(Cd, nbrC);

            forAll(phaseNames_, phasei)
            {
                const label offset = phasei*3*nCells;

                forAll(Cd, celli)
                {
                    displ[offset + 3*celli + d] -= Cd[celli];
                }

                forAll(imageFace_, i)
                {
                    displ[offset + 3*(nInt + i) + d] -= imageCentre_[i][d];
                }

                forAll(coronaFace_, i)
                {
                    displ[offset + 3*(nInt + nImg + i) + d] -=
                        nbrC[coronaFace_[i]];
                }
            }
        }
    }

    rcfd_of_prep_frame_db(seti, vof.begin(), origin.begin(), displ.begin());
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
    coFrame_(),
    nSubFrame_(),
    writeFields_(false),
    labelForm_("auto"),
    displacementForm_(false),
    periodicCentre_(),
    periodicNormal_(),
    periodicShift_(),
    periodicPlaneOfPatch_(),
    periodicImages_(true),
    imageFace_(),
    imageCentre_(),
    periodicDepth_(),
    labels_(),
    extraIntervals_(),
    setInterval_(),
    setFrames_(),
    setDir_(),
    nextFrameTime_(),
    setDone_(),
    opened_(false),
    recording_(false),
    done_(false),
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

    // One label set for the database of the header and one for every
    // further frame interval
    const label nSets = 1 + extraIntervals_.size();

    setInterval_.setSize(nSets, 0);
    setFrames_.setSize(nSets, 0);
    setDir_.setSize(nSets);
    nextFrameTime_.setSize(nSets, great);
    setDone_.setSize(nSets, false);
    coFrame_.setSize(nSets, 0);
    nSubFrame_.setSize(nSets, 1);

    labels_.setSize(3*phaseNames_.size()*nSets);

    for (label seti = 0; seti < nSets; seti++)
    forAll(phaseNames_, phasei)
    {
        for (direction d = 0; d < 3; d++)
        {
            labels_.set
            (
                labelIndex(seti, phasei, d),
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            word
                            (
                                "rCFD" + cmptNames[d]
                              + (seti ? Foam::name(seti) : word::null)
                            ),
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

    setBoundaryFaces();
    setPeriodic();

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

    extraIntervals_ =
        dict.lookupOrDefault<scalarList>("extraFrameIntervals", scalarList());

    if (extraIntervals_.size() > RCFD_OF_MAX_DB - 1)
    {
        FatalIOErrorInFunction(dict)
            << "extraFrameIntervals holds " << extraIntervals_.size()
            << " intervals, at most " << RCFD_OF_MAX_DB - 1 << " are possible"
            << exit(FatalIOError);
    }

    forAll(extraIntervals_, i)
    {
        if (extraIntervals_[i] <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "extraFrameIntervals: " << extraIntervals_[i]
                << " s is not a frame interval" << exit(FatalIOError);
        }
    }

    labelForm_ = dict.lookupOrDefault<word>("labelForm", "auto");
    periodicImages_ = dict.lookupOrDefault<Switch>("periodicImages", true);

    if
    (
        labelForm_ != "auto"
     && labelForm_ != "position"
     && labelForm_ != "displacement"
    )
    {
        FatalIOErrorInFunction(dict)
            << "labelForm must be auto, position or displacement"
            << exit(FatalIOError);
    }

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
    const scalar tol = 1e-6*min(setInterval_);

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

            forAll(nextFrameTime_, seti)
            {
                nextFrameTime_[seti] = tOld + setInterval_[seti];
            }

            Info<< type() << " " << name() << ": recording starts at t = "
                << tOld << " s (start of the run)" << endl;
        }
        else
        {
            if (t >= recordStart_ - tol)
            {
                resetLabels();

                recording_ = true;

                forAll(nextFrameTime_, seti)
                {
                    nextFrameTime_[seti] = t + setInterval_[seti];
                }

                Info<< type() << " " << name() << ": recording starts at t = "
                    << t << " s" << endl;
            }

            return true;
        }
    }

    // Advance the labels of all sets over the time step just solved
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

    forAll(setDone_, seti)
    {
        if (setDone_[seti]) continue;

        const scalar tolSet = 1e-6*setInterval_[seti];

        // "db_10ms/data" for a further database
        const std::string dbName
        (
            seti ? std::string(setDir_[seti]) + "data" : std::string("data")
        );

        if (t < nextFrameTime_[seti] - tolSet) continue;

        if (t > nextFrameTime_[seti] + tolSet)
        {
            nLateFrames_++;

            if (nLateFrames_ <= 5)
            {
                WarningInFunction
                    << "Frame taken at t = " << t << " s, "
                    << t - nextFrameTime_[seti] << " s after its time. Use "
                    << "adjustTimeStep yes or a time step that divides the "
                    << "frame interval"
                    << (seti ? " (" + dbName + ")" : std::string()).c_str()
                    << "." << endl;
            }
        }

        Info<< type() << " " << name() << ": frame "
            << rcfd_of_prep_frames_written_db(seti)
            << (seti ? " of " + dbName : std::string()).c_str()
            << ": label Courant number up to "
            << coFrame_[seti] << ", up to " << nSubFrame_[seti]
            << " label sub-steps per time step"
            << (
                   coFrame_[seti]
                 > labelMaxCo_*maxLabelSubCycles_*(1 + 1e-6)
                 ? " (limited by maxLabelSubCycles)"
                 : ""
               )
            << endl;

        coFrame_[seti] = 0;
        nSubFrame_[seti] = 1;

        handOverFrame(seti);
        resetLabels(seti);

        nextFrameTime_[seti] += setInterval_[seti];

        if (rcfd_of_prep_frames_written_db(seti) >= setFrames_[seti])
        {
            // This database is complete; the others go on
            rcfd_of_prep_close_db(seti);

            setDone_[seti] = true;

            if (setDone_.size() > 1)
            {
                Info<< type() << " " << name() << ": "
                    << dbName.c_str()
                    << " complete at t = " << t << " s ("
                    << setFrames_[seti] << " frames of "
                    << setInterval_[seti] << " s)" << endl;
            }
        }
    }

    bool allDone = true;
    forAll(setDone_, seti)
    {
        allDone = allDone && setDone_[seti];
    }

    if (allDone)
    {
        rcfd_of_prep_close();

        done_ = true;

        Info<< type() << " " << name() << ": database complete at t = "
            << t << " s (" << nFrames_ << " frames";

        if (setDone_.size() > 1)
        {
            Info<< " and " << setDone_.size() - 1 << " further databases";
        }

        Info<< ", label transport "
            << cpuLabels_ << " s cpu on the master). The simulation can "
            << "be stopped." << endl;
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

    // The next frame of any database that is still recording
    scalar tNext = vGreat;

    forAll(setDone_, seti)
    {
        if (!setDone_[seti])
        {
            tNext = min(tNext, nextFrameTime_[seti]);
        }
    }

    return max(tNext - t, scalar(0));
}


// ************************************************************************* //
