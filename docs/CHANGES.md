# Changes

Newest first. One entry per piece of work: date, topic, files, what changed, tested / not tested.

## 2026-11-23  Several granulation runs in one picture (granulation chat)

- `py/rcfd_granulation_compare.py` (new): overlays runs or batches in stages (size, span, variance,
  difference between the chambers) and plots variance added per minute and run time per process
  minute against the frame interval; table as CSV. For the study of the frame interval on the GF25
  batch. Reads the run folders only.
- Tested: four runs on the 10, 20, 25 and 30 ms databases of one small fluidised-bed recording, and
  a single run against a batch in three stages. Not tested: GF25 runs.
- GF25 15 kg recording with `extraFrameIntervals (0.010 0.015 0.020 0.025)` done on the cluster
  (r7): D_C2C of the solid 4.0e-3, 3.1e-3, 2.5e-3, 2.0e-3, 1.7e-3 m2/s at 5 to 25 ms. The solid
  moves less than one cell per frame (0.17 to 0.87 h), so D_C2C falls more slowly than 1/interval.

## 2026-11-20  Release r7: several frame intervals in one recording (granulation chat)

- `functionObject/rCFDPrep.H/.C`: entry `extraFrameIntervals`; one label set per interval
  (`labels_[3*(set*nPhases + phase) + d]`), frame times, completion and messages per set.
- `src/rcfd_of_bridge.h/.c`: one library instance per database (`bridge_db_t`, `bridge_select`,
  `bridge_store`); settings `n_extra`, `extra_intervals`; info `n_db`, `db_frames`, `db_interval`,
  `db_dir`; new calls `rcfd_of_prep_frame_db`, `rcfd_of_prep_frames_written_db`,
  `rcfd_of_prep_close_db`. The old calls act on database 0. A user function `user_set_file_dict` is
  now called through the bridge. `src/rcfd_of_images.h`: paths with the folder of the database.
- `lib/`, `src/rcfd_run_main.c`, `user/`: unchanged.
- `tutorial/GF25_fluidized/run_gf25_batch.sbatch`: frame interval and frames of a stage from its
  database.
- Tested: see manual 2.9, release r7. Not tested: cluster, more than three further intervals.
- Planned next, not started: second propagator moment (six more label fields per phase) for a local
  diffusivity `data/local_diff`, as an option per phase.

## 2026-11-20  Recurrence statistics of a database (granulation chat)

- `py/rcfd_recurrence.py` (new): reads `data/` and writes `recurrence/` in the case folder:
  recurrence matrix and distance against time lag, seam of the replay loop, random-cell series with
  autocorrelation and spectrum, bed height and phase per chamber, plane maps, displacement per frame
  and numerical diffusivity D_C2C of the shifts. No change to `lib/`, `src/` or the run.
- `docs/manual.md` 1.5: section "Recurrence statistics of a database", figure 14.
- Tested: small fluidised bed (40 and 400 frames, 4 ranks), GF25 mesh with a two-frame database.
  Not tested: a full GF25 database (2000 frames, 32 ranks) and meshes with periodic sides (the
  displacement across a periodic side is counted as a long one).

## 2026-11-19  Release r6: runs in stages on several databases (granulation chat)

- `src/rcfd_run_main.c`: option `-init <folder>` (start from `out/state_final.bin` of another run on
  the same mesh, cell by cell through the global cell numbers), state file written at the end of every
  run, stop request from the user header. Snapshots unchanged.
- `user/udf.h`, `user/rcfd_user_wrap.c`: `rcfd_init_folder`, `rcfd_init_vof`, `rcfd_init_time`,
  `rcfd_request_stop`.
- `tutorial/GF25_fluidized/rCFD/rCFD_user_run_granulation.h`: batch granulation with three moment
  fields; particle number from the start mass; stage limit from `rCFD/stage_mass_limit.txt`; totals in
  `rCFD_granulation_state.txt`; cells that held no solid before a switch get the mean over all
  particles; joint correction of the three fields. Defaults `_face_swap_diffusion 0.0` and
  `_spray_height 0.05` are a provisional calibration on one batch.
- `tutorial/GF25_fluidized/run_gf25_batch.sbatch`: chains the stages listed in `stages.txt`.
- `py/rcfd_granulation_plot.py`: several stage folders, variance rate, switches marked.
- Tested: two small databases of one mesh, three stages and six switches without spray (totals kept);
  `fluidisedBed` mixing tutorial gives the same monitor as before. Not tested: GF25 on the cluster.

## 2026-11-17  GF25: batch granulation with three moment fields

- Header, plot script, equation pictures, `docs/notes/moments_vs_bins_test.py`.

## Earlier

See `git log` and manual 2.9 (releases r2 to r5).
