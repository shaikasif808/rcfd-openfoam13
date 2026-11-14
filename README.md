# rCFD inside OpenFOAM-13

Same workflow as with Fluent: two user headers per case, compile, prepare, run.

## This repository

| Folder | Content |
|---|---|
| `lib/` | rCFD library, rCFD-PFM 26.03.0 with the TUHH changes, unchanged for OpenFOAM |
| `src/`, `functionObject/`, `user/`, `bin/` | the layer between OpenFOAM-13 and the library: coupling, function object `rCFDPrep`, run program, user wrapper, `rcfd_compile` |
| `py/` | `snap_to_foam.py`: rCFD snapshots to OpenFOAM time directories for ParaView |
| `tutorial/fluidisedBed` | multiphaseEuler tutorial, size classes with drift |
| `tutorial/aeratedStirredTank` | gas-liquid tank with MRF, species in the gas |
| `tutorial/GF25_fluidized` | GF25 fluidised bed, Euler-Euler with kinetic theory, with geometry and job script |
| `docs/manual.md` | the manual: use, what was built and why, internals, change log |

The text below and the manual call the package folder `rcfd_foam2`. Clone under that name:

```
git clone <url of this repository> rcfd_foam2
```

Private repository. The library in `lib/` is the work of the PFM group at JKU Linz
(S. Pirker) with changes made at TUHH: do not make the repository public without their
agreement.

| Fluent | Here |
|---|---|
| `rCFD_user_prep.h` in the user source folder | `rCFD/rCFD_user_prep.h` in the case |
| `rCFD_user_run.h` | `rCFD/rCFD_user_run.h` in the case |
| build the UDF library | `rcfd_compile` in the case directory |
| `rcfd_prep.scm` (CFD run with the prep UDF) | `foamRun` with the entry `rCFDPrep1` in `system/functions` |
| `rcfd_run.scm` | `mpirun -np N rCFD/rcfd_run -steps 800 -snap 10` |
| contours of user-defined memories | `rcfd_user_field()` in `user_post`, `snap_to_foam.py`, ParaView |

The user headers have the same functions as in Fluent (`user_set_solver_dict`,
`user_set_phase_dict`, `user_set_data_dict`, `user_set_balance_dict`, `user_init_data`,
`user_access_data_before_shift`, `user_access_data_after_swap`, `user_post`, ...,
`user_register_functions`) and the same first lines (`#include "udf.h"`,
`#include "rCFD_library.h"`). The rCFD library in `lib/` is rCFD-PFM 26.03.0 with the TUHH
changes of `src_opt_fable_mpifix` (md5 of `rCFD_library.c`:
6fa34fb96a9fa47d1954fbf622821d07); nothing in it was changed for OpenFOAM.

## Where to edit what

| What | Where |
|---|---|
| frame interval (`global_time_step`), number of frames, number of phases, interpolation, vof cutoff, recurrence format | `rCFD/rCFD_user_prep.h` |
| which OpenFOAM phases are recorded, start of the recording | `system/functions`, entry `rCFDPrep1` (`phases`, `recordStart`) |
| number of species (data) per phase, densities, drift, diffusion, balances, initial state, sources, monitors, output arrays | `rCFD/rCFD_user_run.h` |
| number of rCFD time steps, snapshot interval | command line of `rCFD/rcfd_run` (`-steps`, `-snap`) |
| CFD itself (end time, time step, write interval) | `system/controlDict` as usual |

`global_time_step`, `number_of_frames` and `number_of_phases` appear in both headers, as in
Fluent; the run program stops with a message if they do not fit the database.

## Build (once)

OpenFOAM-13 environment sourced, `mpicc` of the same MPI in the PATH:

```
cd rcfd_foam2
./Allwmake
```

Writes `librcfdOF.so` and `librCFDFunctionObjects.so` to `$FOAM_USER_LIBBIN` and
`rcfd_compile` to `$FOAM_USER_APPBIN`. Do not move the folder afterwards: `rcfd_compile`
uses its sources.

## A case, step by step

```
cd <case>                         # an OpenFOAM case, e.g. a copy of the fluidisedBed tutorial
rcfd_compile -templates           # copies the two example headers into ./rCFD (never overwrites)
```

1. Edit `rCFD/rCFD_user_prep.h` (frame interval, number of frames).
2. Add to `system/functions`:

   ```
   rCFDPrep1
   {
       type            rCFDPrep;
       libs            ("librCFDFunctionObjects.so");

       phases          (air particles);    // order = phase index in the user headers
       recordStart     2.0;                // s, end of the spin-up
   }
   ```

   `endTime` in `controlDict` at least `recordStart` + frames x frame interval.
3. `rcfd_compile prep`, then run the CFD as usual (`foamRun`, or `decomposePar` and
   `mpirun -np N foamRun -parallel`). The database is written to `./data` during the run;
   the log says `database complete`.
4. Edit `rCFD/rCFD_user_run.h` (species, model), `rcfd_compile run`, then

   ```
   mpirun -np N rCFD/rcfd_run -steps 800 -snap 10      # N as in the CFD run
   python3 <path>/rcfd_foam2/py/snap_to_foam.py --run . --case . --out rCFD_results
   ```

   Step 4 can be repeated as often and as long as wanted on the same database: change the
   header, `rcfd_compile run`, run again. No CFD is needed for that.

After every change of a header: `rcfd_compile` (both), `rcfd_compile prep` or
`rcfd_compile run`. `foamRun` (entry `rCFDPrep1`) and `rCFD/rcfd_run` stop with a message
if the header next to the build was changed after it was compiled (cksum of the header,
stored in the build), so a forgotten compile cannot go unnoticed.

## Preparation: notes

* The position labels of each phase are transported every solver time step with the
  phase flux of the solver (`alphaRhoPhi.<phase>`), starting from the cell centres, and set
  back after every frame. With `adjustTimeStep yes` the time step is cut so that the frame
  times are met exactly; with a fixed time step choose one that divides the frame interval.
  No `executeControl`/`writeControl` in the entry.
* `shifts_from_displacement_fields_on` and `time_steps_per_monitoring_interval` are set by
  the coupling; `fluent_set_cell_dict`/`fluent_set_face_dict` in the header take the cell
  and face counts of the OpenFOAM mesh.
* One log line per frame:
  `frame 12/80: phase 0: ... | phase 1: 4553 primary, 13649 total shifts (114 from other ranks), displacement rms/max 5.35/13.68 mm`.
  In a parallel run "from other ranks" must not be zero where material crosses processor
  boundaries.
* A restart of the CFD at `recordStart` (spin-up run, then database run) records from the
  first time step.
* Other solvers: any solver that has a flux per recorded phase. Names per phase, and a
  single-phase solver without phase fraction:

  ```
  phases          (fluid);
  fields
  {
      fluid { alpha none; flux phi; rho none; }
  }
  ```

  Defaults are `alpha.<phase>`, `alphaRhoPhi.<phase>`, `rho.<phase>` (multiphaseEuler).
* Optional entries: `userLibrary` (`rCFD/librCFD_user_prep.so`), `normPhase` (last phase),
  `divScheme` (`"Gauss limitedLinear 1"`), `solver { }`, `residualAlpha` (1e-6),
  `writeFields` (no; yes writes `rCFDdisplacement.<phase>` at write times).
  `labelMaxCo` (1), `maxLabelSubCycles` (200): see the next point.
* Large solver Courant numbers (rotating zones with MRF, jets, `maxCo` far above 1): the
  labels are advanced in sub-steps with the flux of the time step, so that their Courant
  number stays below `labelMaxCo`. The log gives the label Courant number and the number of
  sub-steps per frame. In the vortex test a solver step of a whole frame (Courant 6.5) gives
  7 to 14 % label error without sub-steps, 1 % with `labelMaxCo 0.5` and 0.2 % with 0.1.
  The cost grows with the number of sub-steps; before `recordStart` there is none.
* Inside an MRF zone the solver transports with the flux relative to the rotor, and so do the
  labels: material there moves by the rotor speed relative to the mesh (54 to 75 degrees per
  25 ms at 500 rpm in the stirred tank). Use `labelMaxCo 0.5` or smaller there.
* A mesh made in parallel (`snappyHexMesh -parallel`) has no cell addressing: run
  `reconstructPar` once before `foamRun` (the function object says so if it is missing).
* The boundary patches are written to `data/topo/patches_<rank>` for the run mode
  (`rcfd_patch_faces`, below).
* A run that ends before all frames are recorded (end time, or a change of `controlDict`
  while it runs) leaves a database with the frames written and a replay path; the log says
  which `number_of_frames` to set in `rCFD_user_run.h`.

## Run: notes

* `rcfd_run` holds no model. Everything comes from `rCFD/rCFD_user_run.h`.
* Output: `out/snap_%06d.bin` with the volume fraction of every phase, every data of every
  phase and the arrays made with `rcfd_user_field()`; names in `out/fields.txt`.
  `snap_to_foam.py` writes them as `rCFD_vof.<phase>`, `rCFD_c<i>.<phase>`, `rCFD_<name>`
  into `rCFD_results/<time>/`; open `rCFD_results/rCFD_results.foam` in ParaView.
  ParaView hides time 0 by default: untick "Skip Zero Time" in the Properties of the
  `.foam` file to see the initial state.
  Monitors written by `user_post` (e.g. `rCFD_mean_height_monitor.out`) and the library's
  `rCFD_balance.out` as in Fluent.
* Boundaries in the run header (the counterpart of a face loop over a boundary thread):

  ```
  const int *c; int i, n = rcfd_patch_faces(handle, "inlet", &c, NULL, NULL);
  for(i = 0; i < n; i++) data_gas[c[i]][tracer] = 1.0;      /* fixed value at an inlet */
  ```

  `n` is the number of faces of the patch on this processor (-1: no such patch), `c[i]` the
  cell next to face i; the fourth and fifth argument return face areas and face centres.
  The cells at an inlet take their content from themselves in the shift, so a value set
  there every step (`user_access_data_before_shift`) acts as the inlet condition.
  A species that enters or leaves the domain is not conserved: set the target of its mass
  balance to the current mass every step (`rcfd_set_balance_mass_integral_target` and
  `..._target_global` in `user_access_data_after_swap`), otherwise the balance correction
  pulls the mass towards the start value (in the test about half of the species fed was
  taken out again). `tutorial/aeratedStirredTank` shows both.
* The example header: gas with one data, solids with `_number_of_particle_fractions` size
  classes and the drift law `u_i = -k_d (d_i - d32)|d_i - d32|`. Three classes instead of
  two: `_number_of_particle_fractions 3` and three entries in `Diameter[]` and
  `Mass_fraction[]`. `_tracer_height` > 0 starts a mixing test instead.
  Check the drift velocity before a long run: with 250 / 350 um, d32 300 um and
  k_d = 2e5 it is 0.5 mm/s and nothing segregates in a bubbling bed; the 2e5 of the
  example belongs to its 1.8 / 3.3 mm classes (0.03 and 0.24 m/s).

## Bringing a Fluent user header over

Copy it to `rCFD/rCFD_user_prep.h` or `rCFD/rCFD_user_run.h` and run `rcfd_compile`.
`udf.h` of this package (in `user/`) provides `RP_NODE`, `myid`, `Message`, `Message0`,
`PRF_GRSUM1`, `PRF_GRHIGH1`, `PRF_GRLOW1`, `PRF_GISUM1`, `REGISTER_RCFD_UDF`,
`fluent_set_cell_dict`, `fluent_set_face_dict`. What has no counterpart:

* loops over Fluent threads with `C_UDMI` (and `Get_Domain`, `thread_loop_c`, ...): put
  them between `#ifndef RCFD_OPENFOAM` and `#endif`, and fill a named array instead:

  ```
  double *vf = rcfd_user_field(handle, "vf_diam_0");       /* length: internal cells */
  if(vf) for(i_cell = 0; i_cell < int_cell_count; i_cell++) vf[i_cell] = ...;
  ```

* user functions the library version in `lib/` does not know (for example
  `user_set_calibration_dict`): they are reported and ignored.

## Tutorial

```
cd rcfd_foam2/tutorial/fluidisedBed
./Allrun            # 4 processors;  ./Allrun serial
```

Copies `$FOAM_TUTORIALS/multiphaseEuler/fluidisedBed` into this folder (if there is no
`system` yet), sets `endTime 4` and `writeInterval 0.5`, appends `rCFDPrep.functions` to
`system/functions`, then `rcfd_compile`, `foamRun`, `rCFD/rcfd_run -steps 800 -snap 10`
and `snap_to_foam.py`. The headers are in `rCFD/` of that folder.

## Second tutorial: a species in the gas of an aerated stirred tank

```
cd rcfd_foam2/tutorial
cp -r aeratedStirredTank aeratedStirredTank_test
cd aeratedStirredTank_test && ./Allrun test       # the chain on 8 frames of the start-up
cd ../aeratedStirredTank && ./Allrun              # the tutorial to 15 s, 40 frames from 14 s
```

Copies `$FOAM_TUTORIALS/multiphaseEuler/aeratedStirredTankMRF` (gas and liquid, MRF, mesh by
`snappyHexMesh`), records both phases and transports a species in the gas that is fed at the
patch `inlet` (`rCFD/rCFD_user_run.h`: `rcfd_patch_faces`, open mass balance, monitor
`rCFD_gas_tracer_monitor.out`, field `rCFD_vf_tracer` = species x gas volume fraction).

## Third case: the GF25 fluidised bed

`tutorial/GF25_fluidized` is a complete case (dictionaries, geometry, mesh dictionaries, drag
model, user headers, job script), not a copy of an OpenFOAM tutorial. Steps and the values
that are still assumed are in `tutorial/GF25_fluidized/README.txt`.

* 15 kg of 1.8 mm particles (1040 kg/m3), 2.135 m/s fluidisation air, four nozzles.
* Mesh by `blockMesh` and `snappyHexMesh` from `constant/triSurface`: 177 333 cells.
* Drag `BeetstraLimited` (`dragModel/`, built with `wmake libso`): the Beetstra model of
  OpenFOAM-13 with a lower limit on the Reynolds number, because the original divides by
  zero where both phases have the same velocity.
* 2000 frames of 5 ms from 5 s. `run_gf25.sbatch` runs all steps on 32 processes; the
  folder in its `--chdir` line is that of the TUHH cluster.
* Tested for the first 20 ms on 4 processes (meshing, start-up, frames, `rcfd_run`). The
  fluidised state has not been simulated yet.

## Tested (OpenFOAM-13 built from the 13-master source, Open MPI 4.1, gcc 13)

* `wmake` and `rcfd_compile` without errors or warnings.
* Run program with the example header against the earlier fixed two-fraction program on the
  same database: snapshots identical (difference 0).
* Three size classes by editing the example header; single-phase preparation and run with a
  volumetric flux (`alpha none; flux phi;`).
* Label transport against an exact back-trace in a prescribed vortex: 0.2-0.3 % error of
  the displacement; in the fluidisedBed tutorial against a back-trace through velocity files
  written every 5 ms: 18.42 mm against 18.46 mm rms, correlation 0.98-0.99.
* fluidisedBed tutorial on 4 processors, 80 frames of 0.025 s: a tracer that starts in
  processor 0 only reaches the same mass fraction in all four processors; mass error below
  0.1 % with the balance correction.

* aeratedStirredTankMRF (254 464 cells from snappyHexMesh, MRF at 500 rpm, solver Courant
  number up to 36, 4 processors), 8 frames of 0.025 s during the start-up with gas fed from
  the start of the recording: the cell numbering of the snapshots equals that of the
  reconstructed mesh (gas fraction against OpenFOAM's own field: 8e-8); the gas tagged at
  the inlet grows by 0.43 to 0.47 L per step against 0.467 L fed by the CFD, 79 % of the
  plume cells are above 0.9 after 7 steps, and the head-space gas drawn down at the shaft
  stays untagged. Without label sub-steps the labels in the rotating zone were useless
  (origin points up to 12 m away). Not tested: the steady state after 15 s, a long run.

## Limits

* Static mesh. Cyclic and other coupled patches apart from processor patches are refused.
* The run must use the decomposition (number of processors) of the CFD run that wrote the
  database. Every processor needs at least one cell and one face.
* A restart of the CFD inside the recording window starts the database again.
* Run-mode output goes through snapshots and `snap_to_foam.py`, not through OpenFOAM's own
  field output.
* Large displacements cost time in the preparation: where the origin of a cell is not among
  its neighbours the library searches all cells of the processor (stirred tank: 56 to 97 s
  per frame on 2 cores for 250 000 cells, more than the CFD of that frame).
* Library, not changed here: the preparation keeps some memory per frame
  (`rCFD_parallel.h:209`), about 1 GB for 5e5 cells, 16 ranks and 320 frames.
