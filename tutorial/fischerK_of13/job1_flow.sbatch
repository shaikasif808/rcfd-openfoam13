#!/bin/bash -l
#SBATCH --job-name=fischerK_flow
#SBATCH --chdir=/work/vt1/csa0844/25_rCFD_openfoam_tut/fischerK_of13
#SBATCH --constraint=OS9
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=10000
#SBATCH --time=12:00:00
#SBATCH --output=log.slurm.%j
#
# Job 1 of the 16 million cell lattice case in OpenFOAM-13: mesh and steady flow.
#
#   1. mesh      points, faces, owner, neighbour and boundary are COPIED from the v2512 case
#                (3.6 GB). Nothing is written to the v2512 case.
#   2. sides     the four cyclicAMI patches become ordinary patches and are coupled in pairs
#                with createNonConformalCouples
#   3. decomposePar, checkMesh
#   4. foamRun   steady flow (incompressibleFluid), up to 3000 iterations
#
# Steps that are done are skipped: the script can be submitted again and the flow then
# continues from the last written iteration.
#
# Disk space in this folder: about 10 GB (mesh 3.6 GB, 64 processor folders about 3 GB,
# two flow states of about 1 GB each, logs).

# --- source case and environment -----------------------------------------------
SRC=${SRC:-/work/vt1/csa0844/25_rCFD_openfoam_tut/fischerK_a1_0_b1_0_c1_0}

if [ -z "$RCFD_NO_MODULES" ]; then          # the lines of the of13 function
    conda deactivate 2>/dev/null
    module purge
    module load openmpi/5.0.7
    source "$HOME/OpenFOAM/OpenFOAM-13/etc/bashrc"
fi

stop() { echo "STOP: $*"; exit 1; }

echo "job ${SLURM_JOB_ID:-none} on $(hostname), $(date)"
echo "case $(pwd)"

command -v foamRun > /dev/null                   || stop "OpenFOAM-13 is not loaded"
command -v createNonConformalCouples > /dev/null || stop "createNonConformalCouples not found in this OpenFOAM build"
[ -f system/controlDict ]                        || stop "no system/controlDict in $(pwd)"

# --- 1. mesh -----------------------------------------------------------------------
pm=constant/polyMesh

if [ ! -f $pm/boundary ]; then

    for f in points faces owner neighbour boundary; do
        [ -f "$SRC/$pm/$f" ] || stop "$SRC/$pm/$f is missing"
    done

    echo "copying the mesh from $SRC ($(du -shc "$SRC/$pm"/{points,faces,owner,neighbour} | tail -1 | cut -f1))"

    mkdir -p $pm
    for f in points faces owner neighbour; do
        cp "$SRC/$pm/$f" $pm/ || stop "copy of $f failed"
    done

    # boundary: a cyclicAMI patch becomes an ordinary patch, its AMI entries go
    awk '
        /^[ \t]*inGroups[ \t].*cyclicAMI/                 { next }
        /^[ \t]*type[ \t]+cyclicAMI;/                     { sub(/cyclicAMI/, "patch"); ami = 1; print; next }
        ami && /^[ \t]*}/                                 { ami = 0 }
        ami && /^[ \t]*(matchTolerance|transform|neighbourPatch|separationVector|AMIMethod|restartUncoveredSourceFace|fraction)[ \t]/ { next }
        { print }
    ' "$SRC/$pm/boundary" > $pm/boundary

    grep -q cyclicAMI $pm/boundary && stop "cyclicAMI is still in $pm/boundary"
fi

n_cells=$(grep -m1 -o "nCells: *[0-9]*" $pm/owner | grep -o "[0-9]*$")
echo "mesh: ${n_cells:-?} cells, patches: $(grep -B2 "^ *type " $pm/boundary | grep -v "type\|{\|--" | tr -s ' \n' ' ')"

# --- 2. periodic sides -------------------------------------------------------------
if ! grep -q nonConformalCyclic $pm/boundary; then
    createNonConformalCouples -overwrite > log.createNonConformalCouples 2>&1 \
        || stop "createNonConformalCouples failed, see log.createNonConformalCouples"
fi
echo "coupling of the sides (coverage 1/1/1 = the two sides cover each other completely):"
grep -E "source patch|target patch|coverage" log.createNonConformalCouples | sed 's/^ */    /'

# --- 3. decomposition ----------------------------------------------------------------
if [ ! -f "$FOAM_LIBBIN/libscotchDecomp.so" ] && grep -q "method *scotch" system/decomposeParDict; then
    echo "no scotch library in this build: hierarchical decomposition"
    cp system/decomposeParDict.hierarchical system/decomposeParDict
fi

np=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)
[ "$np" = "${SLURM_NTASKS:-$np}" ] || stop "numberOfSubdomains is $np but the job has $SLURM_NTASKS tasks"

if [ ! -d processor0 ]; then
    rm -rf 0; cp -r 0.orig 0
    decomposePar > log.decomposePar 2>&1 || stop "decomposePar failed, see log.decomposePar"
    grep -E "Max number of cells|Number of processor faces" log.decomposePar | sed 's/^ */    /'
fi

if [ ! -f log.checkMesh ]; then
    mpirun -np "$np" checkMesh -parallel > log.checkMesh 2>&1
fi
grep -q "Mesh stats" log.checkMesh || echo "checkMesh did not run, see log.checkMesh (the flow starts anyway)"
grep -E "^ *cells:|Max aspect ratio|Min volume|Mesh non-orthogonality|Max skewness|Mesh OK|Failed [0-9]* mesh" log.checkMesh | sed 's/^ */    /'

# --- 4. steady flow ---------------------------------------------------------------------
log=log.foamRun.flow
[ -f $log ] && log=log.foamRun.flow.$(date +%Y%m%d_%H%M%S)

mpirun -np "$np" foamRun -parallel > "$log" 2>&1
rc=$?

echo "foamRun exit code $rc, $(date)"
grep "^Time = " "$log" | tail -1
grep "Solving for Ux" "$log" | tail -1
grep "Solving for p"  "$log" | tail -1
grep "ExecutionTime" "$log" | tail -1
grep -E "solution converged" "$log" || echo "not converged yet: submit the job again, it continues from the last written iteration"
f=$(ls -t postProcessing/patchAverage*/*/surfaceFieldValue.dat 2> /dev/null | head -1)
[ -n "$f" ] && { echo "mean p/rho at the inlet [m2/s2], last five iterations:"; tail -5 "$f"; }
f=$(ls -t postProcessing/cellMaxMag*/*/volFieldValue.dat 2> /dev/null | head -1)
[ -n "$f" ] && { echo "largest velocity [m/s], last iteration:"; tail -1 "$f"; }
du -sh . | sed 's/^/disk space of the case: /'

exit $rc
