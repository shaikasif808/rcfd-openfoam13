/*---------------------------------------------------------------------------*\
  rCFD preparation mode inside an OpenFOAM run (OpenFOAM-13)       TUHH 11/2026

Class
    Foam::functionObjects::rCFDPrep

Description
    Counterpart of the Fluent UDF of the rCFD preparation mode. While the
    solver runs, position labels (three scalars per phase) are transported
    with the phase flux of the solver, every time step, starting from the cell
    centres. After every frame interval the label of a cell is the position
    its material came from; volume fraction and labels are handed to the rCFD
    library, which builds the cell-to-cell shifts and writes the database
    (./data/c2c, vof, rec, topo, rcfd_db.meta.toml in the case directory).
    The labels are then set back to the cell centres.

    Runs in serial and in parallel on the decomposition of the OpenFOAM run.
    No intermediate files.

    Label equation of phase i (X = one position component), the same as the
    Fluent user-defined scalar, with the continuity error taken out so that a
    uniform label stays uniform:

        ddt(alpha_i rho_i X) + div(alphaRhoPhi_i X)
      - X (ddt(alpha_i rho_i) + div(alphaRhoPhi_i)) = 0

Usage
    \verbatim
    rCFDPrep1
    {
        type            rCFDPrep;
        libs            ("librCFDFunctionObjects.so");

        phases          (air particles);    // order = phase index in rCFD_user_prep.h
        recordStart     2.0;                // s, end of the spin-up
    }
    \endverbatim

    Frame interval (global_time_step), number of frames, interpolation, vof
    cutoff and all other rCFD settings are in the user header of the case,
    rCFD/rCFD_user_prep.h, as with the Fluent UDF. "rcfd_compile prep" builds
    it into rCFD/librCFD_user_prep.so, which is loaded here.

    Do not add executeControl or writeControl: the labels have to be advanced
    every time step. With adjustTimeStep the time step is cut so that the
    frame times are met exactly.

    Optional entries
    \table
        Property      | Description                              | Default
        userLibrary   | library built from rCFD_user_prep.h      | rCFD/librCFD_user_prep.so
        fields        | field names per phase, see below         | alpha.<phase>, alphaRhoPhi.<phase>, rho.<phase>
        normPhase     | phase whose fraction is the recurrence norm | last phase
        divScheme     | convection scheme of the labels          | "Gauss limitedLinear 1"
        solver        | linear solver dictionary                 | PBiCGStab/DILU
        residualAlpha | keeps the label where a phase is absent  | 1e-6
        labelMaxCo    | Courant number above which a solver time step is divided into sub-steps for the labels | 1
        maxLabelSubCycles | largest number of label sub-steps per time step | 200
        writeFields   | write displacement fields at write times | no
    \endtable

    Other field names, or a single-phase solver (no phase fraction):
    \verbatim
        phases          (fluid);
        fields
        {
            fluid { alpha none; flux phi; rho none; }
        }
    \endverbatim

SourceFiles
    rCFDPrep.C

\*---------------------------------------------------------------------------*/

#ifndef rCFDPrep_functionObject_H
#define rCFDPrep_functionObject_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class rCFDPrep Declaration
\*---------------------------------------------------------------------------*/

class rCFDPrep
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Phase names in the phase order of the database
        wordList phaseNames_;

        //- Field names per phase ("none": no such field)
        wordList alphaNames_;
        wordList fluxNames_;
        wordList rhoNames_;

        //- Start of the recording
        scalar recordStart_;

        //- Library built from the user header of the case
        fileName userLibrary_;

        //- Phase whose volume fraction is the recurrence norm (-1: last)
        label normPhase_;

        //- From the user header, read back from the library
        scalar frameInterval_;
        label nFrames_;

        //- Label transport
        string divScheme_;
        dictionary solverDict_;
        scalar residualAlpha_;

        //- Largest Courant number of the label transport: a solver time step
        //  with a larger one is divided into sub-steps for the labels
        scalar labelMaxCo_;

        //- Largest number of label sub-steps per solver time step
        label maxLabelSubCycles_;

        //- Largest label Courant number and sub-step count of the frame
        scalar coFrame_;
        label nSubFrame_;

        //- Write the displacement of the last frame at write times
        Switch writeFields_;

        //- Form of the labels: "position" (the label is the position the
        //  material came from) or "displacement" (the label is the vector
        //  from the cell centre to that position). "auto": displacement if
        //  the mesh has periodic patches, position otherwise
        word labelForm_;
        bool displacementForm_;

        //- Periodic (cyclic) patches: a point beyond the plane (centre,
        //  outward normal) of a patch is moved by the translation to the
        //  neighbour patch
        List<point> periodicCentre_;
        List<vector> periodicNormal_;
        List<vector> periodicShift_;

        //- Plane of every patch (-1: not periodic)
        labelList periodicPlaneOfPatch_;

        //- Cyclic patches without a shift: a non-conformal interface inside
        //  the domain. The cells on its two sides are neighbours like the
        //  cells on the two sides of a face.
        boolList interfacePatch_;

        //- Image cells: behind every face of a periodic patch the library
        //  gets a copy of the cell on the other side, placed where that cell
        //  would be if the mesh continued. They give the library a neighbour
        //  across the periodic patch for the interpolation of the shifts and
        //  are taken out of the database again by the coupling.
        Switch periodicImages_;

        //- Boundary face (index in the list of all patch faces) through
        //  which the values of image cell i are fetched, its
        //  position and the depth of the image layer behind every plane
        labelList imageFace_;
        vectorField imageCentre_;
        scalarList periodicDepth_;

        //- Position labels: [3*phase + direction]
        PtrList<volScalarField> labels_;

        //- State
        bool opened_;
        bool recording_;
        bool done_;
        scalar nextFrameTime_;
        label nLateFrames_;
        scalar cpuLabels_;

        //- Faces of all patches of the finite volume mesh in one list
        //  (for a non-conformal patch these are its coupled faces, not
        //  faces of the polyhedral mesh): start of every patch, the cell
        //  next to every face and its area vector
        labelList bfStart_;
        labelList bfCell_;
        vectorField bfSf_;

        //- Corona cells (cells of neighbour processors that share a face):
        //  boundary face (index in the list above) through which the
        //  values of corona cell i are fetched
        labelList coronaFace_;


    // Private Member Functions

        //- List the faces of all patches of the finite volume mesh
        void setBoundaryFaces();

        //- For every face of a coupled patch (processor, cyclic,
        //  non-conformal cyclic): the value of the cell on the other side
        void neighbourValues
        (
            const scalarField& cellValues,
            scalarList& nbrValues
        ) const;

        //- Find the periodic patches and choose the form of the labels
        void setPeriodic();

        //- Build the topology of this processor and open the library
        void open();

        //- Set the labels to the cell centres
        void resetLabels();

        //- Advance the labels of one phase over the time step just solved
        void solveLabels(const label phasei);

        //- Hand the frame to the library
        void handOverFrame();

        //- Values of the own cells followed by the corona cells
        void cellAndCoronaValues
        (
            const scalarField& cellValues,
            List<double>& values,
            const label stride,
            const label offset
        ) const;


public:

    //- Runtime type information
    TypeName("rCFDPrep");


    // Constructors

        //- Construct from Time and dictionary
        rCFDPrep
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rCFDPrep(const rCFDPrep&) = delete;


    //- Destructor
    virtual ~rCFDPrep();


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Do not execute at the start of the run
        virtual bool executeAtStart() const
        {
            return false;
        }

        //- Advance the labels; hand over a frame when its time is reached
        virtual bool execute();

        //- Write the displacement fields if requested
        virtual bool write();

        //- Close the database if the run ends before it is complete
        virtual bool end();

        //- Time to the next frame, so that the time step meets it
        virtual scalar timeToNextAction();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rCFDPrep&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
