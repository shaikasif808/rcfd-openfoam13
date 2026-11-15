#!/bin/bash -l
#SBATCH --job-name=GF25_rCFD
#SBATCH --chdir=/work/vt1/csa0844/25_rCFD_openfoam_tut/GF25_fluidizedbed_of13
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --constraint=OS9
#SBATCH --time=12:00:00
#SBATCH --output=log.slurm.rcfd.%j
#
# rCFD run of the GF25 fluidised bed on the database in ./data (no CFD).
# Partition, account and memory are not set here: add the lines of your other job scripts
# if the cluster needs them.
#
# Number of rCFD steps and snapshot interval (one step = one frame interval, 5 ms here):
#   sbatch run_gf25_rcfd.sbatch                                   12000 steps = 60 s, snapshot every 1 s
#   sbatch --export=ALL,STEPS=24000,SNAP=400 run_gf25_rcfd.sbatch
STEPS=${STEPS:-12000}
SNAP=${SNAP:-200}
#
# Results of an earlier rCFD run in this folder are not overwritten: they are MOVED to
# ./rcfd_run_<date>_<time>/ before the new run starts (out/, rCFD_results/, log.rcfd_run,
# log.snap_to_foam, rCFD_*.out, rCFD_c2c_step*.txt, rCFD_monitor.png, rCFD_run.trn).
# The database ./data, the CFD results and the headers in ./rCFD are not touched.

# --- environment: the lines of the of13 function ---------------------------------
if [ -z "$RCFD_NO_MODULES" ]; then
    conda deactivate 2>/dev/null
    module purge
    module load openmpi/5.0.7
    source "$HOME/OpenFOAM/OpenFOAM-13/etc/bashrc"
fi

stop() { echo "STOP: $*"; exit 1; }

echo "job ${SLURM_JOB_ID:-none} on $(hostname), $(date)"
echo "case $(pwd)"

command -v rcfd_compile > /dev/null || stop "rcfd_compile not found: run ./Allwmake in the rcfd_foam2 folder"
[ -f data/rcfd_db.meta.toml ]       || stop "no database in ./data (data/rcfd_db.meta.toml is missing)"
[ -f rCFD/rCFD_user_run.h ]         || stop "rCFD/rCFD_user_run.h is missing"

# the run needs the number of processes of the CFD run that wrote the database
np=$(ls data/topo/gid_*.bin 2> /dev/null | wc -l)
[ "$np" -gt 0 ] || stop "no data/topo/gid_*.bin: the database is not complete"
[ "$np" = "${SLURM_NTASKS:-$np}" ] || stop "the database was written by $np processes but the job has $SLURM_NTASKS tasks"

pkg=$(grep -m1 '^RCFD_FOAM_HOME=' "$(command -v rcfd_compile)" | cut -d'"' -f2)

# --- earlier results aside -----------------------------------------------------------
old=""
for f in out rCFD_results log.rcfd_run log.snap_to_foam rCFD_*.out rCFD_c2c_step*.txt rCFD_monitor.png rCFD_run.trn; do
    [ -e "$f" ] && old="$old $f"
done
if [ -n "$old" ]; then
    keep=rcfd_run_$(date +%Y%m%d_%H%M%S)
    mkdir "$keep" && mv $old "$keep"/ || stop "could not move the earlier results to $keep"
    echo "earlier results moved to $keep/:$old"
fi

# --- compile the user headers ----------------------------------------------------------
rcfd_compile > log.rcfd_compile 2>&1 || stop "rcfd_compile failed, see log.rcfd_compile"

# --- pictures while the run is going (needs python3 with matplotlib) -------------------
mon=""
if [ -f "$pkg/py/rcfd_monitor.py" ] && python3 -c "import matplotlib, numpy" 2> /dev/null; then
    python3 "$pkg/py/rcfd_monitor.py" --every 120 > log.rcfd_monitor 2>&1 &
    mon=$!
    echo "rCFD_monitor.png is written every 2 minutes"
else
    echo "no pictures from the job (rcfd_monitor.py or matplotlib not found on this node)."
    echo "On the login node, in this folder:  python3 $pkg/py/rcfd_monitor.py --every 60"
fi

# --- run ---------------------------------------------------------------------------------
echo "rcfd_run: $np processes, $STEPS steps, snapshot every $SNAP steps"
start=$(date +%s)

mpirun -np "$np" rCFD/rcfd_run -steps "$STEPS" -snap "$SNAP" > log.rcfd_run 2>&1
rc=$?

echo "rcfd_run exit code $rc after $(( $(date +%s) - start )) s, $(date)"
grep "run finished" log.rcfd_run
tail -1 rCFD_mean_height_monitor.out 2> /dev/null

[ -n "$mon" ] && kill "$mon" 2> /dev/null
[ -f "$pkg/py/rcfd_monitor.py" ] && python3 "$pkg/py/rcfd_monitor.py" 2> /dev/null

# --- snapshots for ParaView ----------------------------------------------------------------
if [ $rc -eq 0 ] && [ "$SNAP" -gt 0 ]; then
    python3 "$pkg/py/snap_to_foam.py" --run . --case . --out rCFD_results > log.snap_to_foam 2>&1 \
        && echo "ParaView: rCFD_results/rCFD_results.foam" \
        || echo "snap_to_foam.py failed, see log.snap_to_foam (it can be run again on the login node)"
fi

exit $rc
