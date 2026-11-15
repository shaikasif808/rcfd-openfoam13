#!/bin/bash -l
#SBATCH --job-name=fischerK_rCFD
#SBATCH --chdir=/work/vt1/csa0844/25_rCFD_openfoam_tut/fischerK_of13
#SBATCH --constraint=OS9
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=10000
#SBATCH --time=12:00:00
#SBATCH --output=log.slurm.rcfd.%j
#
# Job 2 of the 16 million cell lattice case: rCFD database on the frozen flow of job 1, then
# the rCFD run with diffusion and reactions. Needs package r5 of rcfd_foam2 (./Allwmake).
#
#   1. recording   the converged flow of job 1 becomes time 0 (copied, about 1 GB),
#                  system/ is replaced by system.record/ (the flow settings are kept in
#                  system.flow/), foamRun with the function object rCFDPrep writes two frames
#                  of 0.01 s to ./data. Skipped if the database is complete.
#   2. rCFD run    rCFD/rcfd_run with rCFD/rCFD_user_run.h
#   3. snapshots   snap_to_foam.py for ParaView (rCFD_results/)
#
# Number of rCFD steps and snapshot interval (one step = 0.01 s):
#   sbatch job2_rcfd.sbatch                                   500 steps = 5 s, snapshots at 0 and 5 s
#   sbatch --export=ALL,STEPS=20,SNAP=0 job2_rcfd.sbatch      short test: 20 steps, no snapshot
#   sbatch --export=ALL,STEPS=0 job2_rcfd.sbatch              recording only
STEPS=${STEPS:-500}
SNAP=${SNAP:-$STEPS}
#
# Results of an earlier rCFD run in this folder are MOVED to ./rcfd_run_<date>_<time>/ before a
# new run starts. The database ./data is never removed by this script.
#
# Disk space added by this job: about 9 GB (flow as time 0: 1 GB, database 3 to 4 GB, two
# snapshots 0.8 GB, rCFD_results 2 GB, shift diagnostics rCFD_c2c_step*.txt about 1 GB, which
# can be deleted after the run). Each further snapshot: 0.4 GB + 1 GB in rCFD_results.

if [ -z "$RCFD_NO_MODULES" ]; then          # the lines of the of13 function
    conda deactivate 2>/dev/null
    module purge
    module load openmpi/5.0.7
    source "$HOME/OpenFOAM/OpenFOAM-13/etc/bashrc"
fi

stop() { echo "STOP: $*"; exit 1; }

echo "job ${SLURM_JOB_ID:-none} on $(hostname), $(date)"
echo "case $(pwd)"

command -v foamRun > /dev/null      || stop "OpenFOAM-13 is not loaded"
command -v rcfd_compile > /dev/null || stop "rcfd_compile not found: run ./Allwmake in the rcfd_foam2 folder"

pkg=$(grep -m1 '^RCFD_FOAM_HOME=' "$(command -v rcfd_compile)" | cut -d'"' -f2)
grep -q rcfd_update_corona_field "$pkg/user/udf.h" 2> /dev/null \
    || stop "the rCFD package in $pkg is older than r5: unpack rcfd_foam2_r5.tar.gz and run ./Allwmake"

[ -d processor0 ]                   || stop "no processor0: run job 1 first"
[ -f rCFD/rCFD_user_run.h ] && [ -d system.record ] || stop "rCFD/ or system.record/ is missing: unpack fischerK_of13_job2.tar.gz"

np=$(ls -d processor[0-9]* | wc -l)
[ "$np" = "${SLURM_NTASKS:-$np}" ] || stop "the case is decomposed for $np processes but the job has $SLURM_NTASKS tasks"

# --- 1. recording ----------------------------------------------------------------------
if [ -f data/rcfd_db.meta.toml ] && grep -q "database complete" log.foamRun.record 2> /dev/null; then

    echo "database in ./data is complete: recording skipped"

else

    [ -e data ] && stop "./data exists but is not complete: look at log.foamRun.record, then move or remove ./data yourself"

    if [ ! -f system.flow/controlDict ]; then

        # the converged flow becomes time 0 of the recording
        last=$(ls -d processor0/[1-9]* 2> /dev/null | sed 's|processor0/||' | sort -n | tail -1)
        [ -n "$last" ] && [ -f "processor0/$last/phi" ] || stop "no flow result with phi in processor0 (job 1 not finished?)"
        grep -q "solution converged" log.foamRun.flow* 2> /dev/null || echo "WARNING: no 'solution converged' in log.foamRun.flow*: the flow of iteration $last is used as it is"

        echo "flow of iteration $last is copied to time 0"
        for p in processor[0-9]*; do
            rm -rf "$p/0"; cp -r "$p/$last" "$p/0" || stop "copy in $p failed"; rm -rf "$p/0/uniform"
        done

        mv system system.flow || stop "could not rename system to system.flow"
        cp -r system.record system
        cp system.flow/decomposeParDict system/
    fi

    rcfd_compile > log.rcfd_compile 2>&1 || stop "rcfd_compile failed, see log.rcfd_compile"

    start=$(date +%s)
    mpirun -np "$np" foamRun -parallel > log.foamRun.record 2>&1
    rc=$?
    echo "recording: foamRun exit code $rc after $(( $(date +%s) - start )) s"

    grep -E "periodic patches|image cells behind|labels in" log.foamRun.record | head -4 | cut -c1-200
    grep -E "frame +[0-9]+/[0-9]+:" log.foamRun.record | cut -c1-230
    grep -E "database written|peak resident memory|database complete" log.foamRun.record | cut -c1-230

    grep -q "database complete" log.foamRun.record || stop "the database is not complete, see log.foamRun.record"
fi

du -sh data | sed 's/^/database: /'

[ "$STEPS" -gt 0 ] || { echo "STEPS=0: no rCFD run"; exit 0; }

# --- 2. rCFD run --------------------------------------------------------------------------
old=""
for f in out rCFD_results log.rcfd_run log.snap_to_foam rCFD_*.out rCFD_run.trn; do
    [ -e "$f" ] && old="$old $f"
done
if [ -n "$old" ]; then
    keep=rcfd_run_$(date +%Y%m%d_%H%M%S)
    mkdir "$keep" && mv $old "$keep"/ || stop "could not move the earlier results to $keep"
    echo "earlier results moved to $keep/:$old"
fi

rcfd_compile > log.rcfd_compile 2>&1 || stop "rcfd_compile failed, see log.rcfd_compile"

echo "rcfd_run: $np processes, $STEPS steps, snapshot every $SNAP steps"
start=$(date +%s)

mpirun -np "$np" rCFD/rcfd_run -steps "$STEPS" -snap "$SNAP" > log.rcfd_run 2>&1
rc=$?

echo "rcfd_run exit code $rc after $(( $(date +%s) - start )) s, $(date)"
grep "\[user\]" log.rcfd_run | head -5
grep "\[user\] step" log.rcfd_run | tail -1
grep "run finished" log.rcfd_run | cut -c1-230
echo "outlet (time, A B C D E, then the content of the domain):"
tail -1 rCFD_lattice_monitor.out 2> /dev/null

[ $rc -eq 0 ] || exit $rc

# --- 3. snapshots for ParaView ---------------------------------------------------------------
if [ "$SNAP" -gt 0 ]; then
    python3 "$pkg/py/snap_to_foam.py" --run . --case . --out rCFD_results > log.snap_to_foam 2>&1 \
        && echo "ParaView: rCFD_results/rCFD_results.foam" \
        || echo "snap_to_foam.py failed, see log.snap_to_foam (it can be run again on the login node)"
fi

du -sh . | sed 's/^/disk space of the case: /'
exit 0
