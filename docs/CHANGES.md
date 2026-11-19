# Changes

Newest first. One entry per piece of work: date, topic, files, what changed, tested / not tested.

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
