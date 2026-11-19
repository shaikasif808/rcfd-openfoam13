#!/bin/bash -l
#SBATCH --job-name=GF25_batch
#SBATCH --chdir=/work/vt1/csa0844/25_rCFD_openfoam_tut/GF25_batch_30kgh
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --constraint=OS9
#SBATCH --time=12:00:00
#SBATCH --output=log.slurm.batch.%j
#
# Batch granulation in stages: every stage replays the database recorded at its bed mass, the next
# stage starts from the end of the one before (rcfd_run -init). Needs package r6 of rcfd_foam2.
#
# The folder of --chdir (create it yourself) holds two files:
#   rCFD_user_run.h   the granulation header with the settings of the batch (spray rate, seeds, ...)
#   stages.txt        one line per stage:   <case folder of the database>   <bed mass at which the stage ends [kg]>
#                     e.g.   /work/.../GF25_fluidizedbed_of13   17.5
#                            /work/.../GF25_20kg_of13           22.5
#                            /work/.../GF25_25kg_of13           27.5
#                            /work/.../GF25_30kg_of13           31.0      <- last line: bed mass at the end of the batch
# The script makes stage_1, stage_2, ... there. A stage folder holds links to data/ and constant/ of its
# database (nothing is written into the database folders), a copy of the header, and the results of the
# stage. A finished stage is skipped when the job is submitted again; an unfinished stage folder stops the
# job (look at its log.rcfd_run, then move or remove that folder yourself).
#
#   sbatch run_gf25_batch.sbatch
#   sbatch --export=ALL,STEPS=400000,SNAP=4000 run_gf25_batch.sbatch
STEPS=${STEPS:-400000}      # most steps of one stage; the header ends a stage when its bed mass is reached
SNAP=${SNAP:-4000}          # snapshot every so many steps (0: first and last only)

if [ -z "$RCFD_NO_MODULES" ]; then          # the lines of the of13 function
    conda deactivate 2>/dev/null
    module purge
    module load openmpi/5.0.7
    source "$HOME/OpenFOAM/OpenFOAM-13/etc/bashrc"
fi

stop() { echo "STOP: $*"; exit 1; }

echo "job ${SLURM_JOB_ID:-none} on $(hostname), $(date)"
echo "batch folder $(pwd)"

command -v rcfd_compile > /dev/null || stop "rcfd_compile not found: run ./Allwmake in the rcfd_foam2 folder"

pkg=$(grep -m1 '^RCFD_FOAM_HOME=' "$(command -v rcfd_compile)" | cut -d'"' -f2)
grep -q rcfd_request_stop "$pkg/user/udf.h" 2> /dev/null \
    || stop "the rCFD package in $pkg is older than r6: unpack rcfd_foam2_r6.tar.gz and run ./Allwmake"

[ -f rCFD_user_run.h ] || stop "rCFD_user_run.h (the granulation header) is missing in $(pwd)"
[ -f stages.txt ]      || stop "stages.txt is missing in $(pwd)"
grep -q rcfd_init_folder rCFD_user_run.h || stop "rCFD_user_run.h is not the granulation header with stages"

batch=$(pwd)
k=0; prev=""; folders=""

while read -r db limit rest; do

    case "$db" in ""|\#*) continue ;; esac

    k=$((k + 1)); s="stage_$k"

    [ -f "$db/data/rcfd_db.meta.toml" ] || stop "stage $k: no database in $db/data"
    echo "$limit" | grep -Eq '^[0-9]+([.][0-9]+)?$' || stop "stage $k: '$limit' is not a bed mass in kg"

    np=$(ls "$db"/data/topo/gid_*.bin 2> /dev/null | wc -l)
    [ "$np" = "${SLURM_NTASKS:-$np}" ] || stop "stage $k: the database was written by $np processes but the job has $SLURM_NTASKS tasks"

    if [ -f "$s/out/state_final.bin" ] && grep -q "the user header ended the run" "$s/log.rcfd_run" 2> /dev/null; then

        echo "stage $k ($db, to $limit kg): finished earlier, skipped"

    else

        [ -e "$s" ] && stop "$s exists but is not finished: look at $s/log.rcfd_run, then move or remove $s yourself"

        mkdir -p "$s/rCFD" || stop "cannot create $s"
        ln -s "$db/data" "$s/data"
        [ -d "$db/constant" ] && ln -s "$db/constant" "$s/constant"
        cp rCFD_user_run.h "$s/rCFD/rCFD_user_run.h"
        echo "$limit" > "$s/rCFD/stage_mass_limit.txt"

        cd "$s" || stop "cannot enter $s"

        rcfd_compile run < /dev/null > log.rcfd_compile 2>&1 || stop "stage $k: rcfd_compile failed, see $s/log.rcfd_compile"

        init=""; [ -n "$prev" ] && init="-init ../$prev"

        echo "stage $k: database $db, to $limit kg, $np processes $init"
        start=$(date +%s)

        mpirun -np "$np" rCFD/rcfd_run -steps "$STEPS" -snap "$SNAP" $init < /dev/null > log.rcfd_run 2>&1
        rc=$?

        echo "  rcfd_run exit code $rc after $(( $(date +%s) - start )) s"
        grep "\[user\] switch of the database\|\[user\] bed mass\|\[user\] ERROR\|run finished" log.rcfd_run | cut -c1-200 | sed 's/^/  /'

        cd "$batch" || exit 1

        [ $rc -eq 0 ] || stop "stage $k: rcfd_run failed, see $s/log.rcfd_run"
        grep -q "the user header ended the run" "$s/log.rcfd_run" \
            || stop "stage $k did not reach $limit kg within $STEPS steps: submit again with a larger STEPS after moving $s away"
    fi

    prev="$s"; folders="$folders $s"

done < stages.txt

[ $k -gt 0 ] || stop "stages.txt holds no stage"

echo "all $k stages done, $(date)"

# picture of the whole batch (needs python3 with numpy and matplotlib)
acc=$(grep -m1 '^#define _growth_acceleration' rCFD_user_run.h | awk '{print $3}')
meas=""; [ -f samples.csv ] && meas="--measured samples.csv"
python3 "$pkg/py/rcfd_granulation_plot.py" --case $folders --acceleration "${acc:-1}" $meas --out rCFD_granulation.png \
    || echo "no picture from the job; on the login node:  python3 $pkg/py/rcfd_granulation_plot.py --case $folders --acceleration ${acc:-1}"

du -sh . | sed 's/^/disk space of the batch folder: /'
exit 0
