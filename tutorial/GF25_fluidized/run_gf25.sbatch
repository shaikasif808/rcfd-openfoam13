#!/bin/bash -l
#SBATCH --job-name=GF25_EE
#SBATCH --chdir=/work/vt1/csa0844/25_rCFD_openfoam_tut/GF25_fluidizedbed_of13
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --constraint=OS9
#SBATCH --time=24:00:00
#SBATCH --output=log.slurm.%j
#
# GF25 fluidised bed, multiphaseEuler (OpenFOAM-13) with the rCFD recording.
# Partition, account and memory are not set here: add the lines of your
# stirred-tank job script if your cluster needs them.
#
# Steps that are already done are skipped, so the script can be submitted again:
#   drag library, mesh, rcfd_compile, decomposition and bed, then foamRun.
# A second submission continues from the last written time (startFrom latestTime).
# Continuing in the middle of the recording (after 5 s) has not been tested.

# --- environment: the lines of the of13 function -----------------------------
conda deactivate 2>/dev/null
module purge
module load openmpi/5.0.7
source "$HOME/OpenFOAM/OpenFOAM-13/etc/bashrc"

# The PSM3/libfabric messages on the g nodes are warnings (the stirred-tank job
# ran with them). Not tested: these two lines should silence them on one node.
#export OMPI_MCA_pml=ob1
#export OMPI_MCA_btl=self,sm

stop() { echo "STOP: $*"; exit 1; }

echo "job ${SLURM_JOB_ID:-none} on $(hostname), $(date)"
echo "case $(pwd)"

command -v foamRun > /dev/null      || stop "OpenFOAM-13 is not loaded"
command -v rcfd_compile > /dev/null || stop "rcfd_compile not found: run ./Allwmake in the rcfd_foam2 folder"
[ -f system/controlDict ]           || stop "no system/controlDict in $(pwd)"

np=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)
[ "$np" = "${SLURM_NTASKS:-$np}" ] \
    || stop "numberOfSubdomains is $np but the job has $SLURM_NTASKS tasks"

# --- 1. drag model BeetstraLimited --------------------------------------------
if [ ! -f "$FOAM_USER_LIBBIN/libGF25Drag.so" ]; then
    (cd dragModel && wmake libso) > log.wmake 2>&1
    [ -f "$FOAM_USER_LIBBIN/libGF25Drag.so" ] || stop "drag library not built, see log.wmake"
fi

# --- 2. mesh -------------------------------------------------------------------
if [ ! -f constant/polyMesh/owner ]; then
    [ -f constant/geometry/casing.stl ] || [ -f constant/triSurface/casing.stl ] \
        || stop "geometry is missing (constant/geometry)"
    blockMesh      > log.blockMesh 2>&1      || stop "blockMesh failed"
    snappyHexMesh  > log.snappyHexMesh 2>&1  || stop "snappyHexMesh failed"
fi
grep "Snapped mesh" log.snappyHexMesh 2> /dev/null

# --- 3. rCFD user headers ------------------------------------------------------
rcfd_compile > log.rcfd_compile 2>&1 || stop "rcfd_compile failed, see log.rcfd_compile"

# --- 4. decomposition and bed --------------------------------------------------
if [ ! -d processor0 ]; then
    decomposePar > log.decomposePar 2>&1 || stop "decomposePar failed"
    mpirun -np "$np" setFields -parallel > log.setFields 2>&1 || stop "setFields failed"
fi

# --- 5. simulation with recording ----------------------------------------------
log=log.foamRun
[ -f $log ] && log=log.foamRun.$(date +%Y%m%d_%H%M%S)

mpirun -np "$np" foamRun -parallel > "$log" 2>&1
rc=$?

echo "foamRun exit code $rc, $(date)"
grep "^Time = " "$log" | tail -1
grep "rCFDPrep rCFDPrep1: database complete" "$log" || echo "database not complete"
f=$(ls -t postProcessing/volIntegrate*/*/volFieldValue.dat 2> /dev/null | head -1)
[ -n "$f" ] && tail -1 "$f"

exit $rc
