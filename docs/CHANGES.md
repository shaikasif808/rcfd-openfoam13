# Changes

Newest first. One entry per piece of work: date, topic, files, what changed, tested / not tested.

## 2026-11-20  Recurrence statistics of a database (granulation chat)

- `py/rcfd_recurrence.py` (new): reads `data/` and writes `recurrence/` in the case folder:
  recurrence matrix and distance against time lag, seam of the replay loop, random-cell series with
  autocorrelation and spectrum, bed height and phase per chamber, plane maps, displacement per frame
  and numerical diffusivity D_C2C of the shifts. No change to `lib/`, `src/` or the run.
- `docs/manual.md` 1.5: section "Recurrence statistics of a database", figure 14.
- Tested: small fluidised bed (40 and 400 frames, 4 ranks), GF25 mesh with a two-frame database.
  Not tested: a full GF25 database (2000 frames, 32 ranks) and meshes with periodic sides (the
  displacement across a periodic side is counted as a long one).
- Planned next, not started: several frame intervals in one recording (`rCFDPrep`,
  `src/rcfd_of_bridge.c`) and the second propagator moment for a local diffusivity.

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
