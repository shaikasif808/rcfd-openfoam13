# Hand-over: working on rcfd-openfoam13 in a new chat

Read this first, then `docs/manual.md` (Part 3 explains how rCFD works, Part 1 how the package is
used, Part 2 what was built and why). State of this note: 20 November 2026, package release r7.

## What this repository is

rCFD is the recurrence CFD method of Pirker and Lichtenegger (JKU Linz): a short, expensive flow
simulation is recorded as a database of frames, and passive quantities are then moved through the
recorded flow for long process times at low cost. The C library `rCFD-PFM 26.03.0` was written for
ANSYS Fluent. This repository couples the unchanged library (`lib/`) to OpenFOAM-13.

| Folder | Content |
| --- | --- |
| `lib/` | the rCFD library (C, header-only modules compiled through `rCFD_library.c`), with a few TUHH fixes listed in manual 2.9 |
| `functionObject/` | `rCFDPrep` (C++): runs inside `foamRun`, transports position labels with the solver's flux and hands cell-to-cell shifts to the library, which writes the database `data/` |
| `src/` | coupling: `rcfd_of_bridge.c` (library side of the function object), `rcfd_run_main.c` (the run program), `rcfd_of_common.h`, `rcfd_of_images.h` (periodic sides) |
| `user/` | `udf.h` (what a user header may call, Fluent names kept), `rcfd_user_wrap.c` |
| `bin/rcfd_compile` | builds the two user headers of a case: `rCFD/rCFD_user_prep.h` (recording) and `rCFD/rCFD_user_run.h` (run) |
| `py/` | post-processing: `snap_to_foam.py`, `rcfd_slices.py`, `rcfd_monitor.py`, `rcfd_granulation_plot.py`, `rcfd_recurrence.py`, `of_reader.py` |
| `tutorial/` | `fluidisedBed`, `aeratedStirredTank`, `latticeReactions`, `fischerK_of13` (16 million cells), `GF25_fluidized` (granulator, incl. batch granulation in stages) |
| `docs/` | `manual.md` with figures, this note, `CHANGES.md` |

## How rCFD works here, in five lines

1. Recording (`foamRun` with `rCFDPrep`): per phase three label fields hold the cell-centre position;
   they are transported with the flux for one frame interval; where a cell's label points to is its
   source. The library turns this into weighted cell-to-cell shifts per frame and writes `data/`.
2. Database `data/`: `c2c/` shifts per frame, phase and rank; `vof/` volume fractions per frame;
   `rec/jump_0_0_0` the order of frames; `topo/` cells, faces, neighbours and global cell numbers
   per rank; `rcfd_db.meta.toml`, `rcfd_foam_prep.txt`. Formats: manual 1.3 and 3.3.
3. Run (`mpirun -np N rCFD/rcfd_run -steps .. -snap ..`): no flow is solved. Per step and phase:
   user hook, shift (`c_new(i) = sum_j w_ij c_old(j) / sum_j w_ij`), face-swap diffusion, optional
   drift, global balance correction, user hooks, `user_post`. Manual 3.4.
4. A model lives in the run header as C functions with the Fluent names (`user_set_*_dict`,
   `user_init_data`, `user_access_data_before_phase_loop`, `user_access_data_before_shift`,
   `user_access_data_after_swap`, `user_post`). Extra calls of this coupling are in `user/udf.h`:
   `rcfd_user_field`, `rcfd_patch_faces`, `rcfd_face_distance`, `rcfd_update_corona`,
   `rcfd_update_corona_field`, `rcfd_init_folder`, `rcfd_init_vof`, `rcfd_request_stop`.
5. Output: `out/snap_%06d.bin` (int64 cells, int32 fields, float64 time, then float32 arrays in
   global cell numbering; names in `out/fields.txt`), `out/state_final.bin` (float64, for a later
   run with `-init`), and whatever the header writes.

A database does not have to come from OpenFOAM: anything that writes the files of `data/` in the
formats of manual 3.3 can be replayed by `rcfd_run`. For data that are not a flux field (for example
particle trajectories) the shifts can be built directly: source cell and target cell of every
particle over one frame interval, with weights. The topology files need a mesh (cells, faces,
neighbours, volumes, centres).

## Building and testing

- `./Allwmake` needs the OpenFOAM-13 environment and `mpicc`: it builds `librcfdOF.so`, the run
  program object, `rcfd_compile` and the function object. Only the last step (`wmake libso
  functionObject`) needs OpenFOAM itself; library, coupling and run program are plain C with MPI.
- A new workspace has no OpenFOAM. Building OpenFOAM-13 from source takes long; for work on the run
  side alone, compile the C parts by hand as `Allwmake` and `bin/rcfd_compile` do and test on a small
  database. A small test database can be recorded with `tutorial/fluidisedBed` (6000 cells, about one
  minute on 4 processes) where OpenFOAM is available.
- Tests that must keep passing after a change of `lib/` or `src/`: the `fluidisedBed` tutorial gives
  the same `rCFD_mean_height_monitor.out`; `latticeReactions` gives the outlet values of its README.

## Working rules of the owner (Asif Shaik, TUHH)

- Concise and direct, British English. Ask when unclear. Propose before anything that costs HPC
  time. Do not deviate from the task.
- He runs all SLURM jobs himself and pastes results back. Commands as copy-paste blocks.
- He copies files to the cluster by hand: name the file, where it goes, and give a `grep -c` check
  and an md5. After writing a file into one of his folders, verify that it landed.
- Never delete files outside the scope of the task; ask before anything that affects his storage.
- His e-mail address and GitHub token are private: never put them into files, commits or requests.
- Commits: author `Claude <noreply@anthropic.com>`, messages end with the attribution lines the
  session gives. Chats usually cannot push: hand him a git bundle
  (`git bundle create <file> <last pushed commit>..main`), he pulls it into his clone and pushes.
- Cluster: OpenFOAM-13 with `openmpi/5.0.7`, loaded by the shell function `of13`; package in
  `/fibus/fs0/15/csa0844/OpenFOAM/rcfd_foam2`; cases under `/work/vt1/csa0844/25_rCFD_openfoam_tut`.

## Two chats on one repository

Chats do not see each other. The repository is the only shared state.

- Start of a session: `git log --oneline -15` and read `docs/CHANGES.md` from the top.
- End of a piece of work: add an entry at the top of `docs/CHANGES.md` (date, topic, files, what
  changed, what was tested and what was not), commit, hand over a bundle.
- Keep to your own files where you can. If you must change a shared file (`src/rcfd_run_main.c`,
  `user/udf.h`, `user/rcfd_user_wrap.c`, `lib/`), say so in `docs/CHANGES.md` and keep the change
  small, so the other chat can merge.
- Files in active use by the granulation work (November 2026): `src/rcfd_run_main.c`, `user/udf.h`,
  `user/rcfd_user_wrap.c`, `tutorial/GF25_fluidized/`, `py/rcfd_granulation_plot.py`.

## Known limits (details in manual 2.8)

- The shift is not exactly conservative; the library corrects global balances per field. Fields
  that depend on each other (moments) must be corrected together in the header.
- Interpolation between a source cell and its face neighbours only; short frame intervals over-mix.
- Source search of the recording is brute force per rank (minutes per frame for 250 000 cells per
  rank). Recording memory about 75 bytes per cell of the whole mesh on every rank.
- One database per run; several databases of the same mesh are chained with `-init` (r6). One recording can
  write several databases with different frame intervals (`extraFrameIntervals`, r7).
- Drift and local diffusion across periodic sides are untested.
