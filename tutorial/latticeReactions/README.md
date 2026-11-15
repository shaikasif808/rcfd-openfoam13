# latticeReactions

Question of the tutorial: if the flow field is known, how well does rCFD reproduce transport,
diffusion and reactions on it, and what does it cost?

## Case

Liquid (density 992.25 kg/m3, viscosity 0.6527 mPa s) flows at 1 mL/s through a Fischer-Koch S
lattice in a channel of 10 x 10 mm, periodic in x and y. Species A enters at the inlet.

| Reaction | Rate |
|---|---|
| A -> B at the lattice wall (patch `reactingWall`) | 1e6 1/s in the cells at the wall: instantaneous, limited by the transport to the wall |
| B -> C | 10 1/s |
| C -> D | 1 1/s |
| 2 B -> E | 55.125 1/s x Y_B^2 |

All species diffuse with D = 6.578e-7 m2/s (Schmidt number 1, as the reacting solver of the
original case does with `muEff`). The mesh is a stair-step mesh of 209 104 cubes of 0.2 mm,
written by `makeLatticeMesh.py`; it stands for a body-fitted mesh, the method is the same.

## Run

```
./Allrun                 # all steps, 4 processors, about 8 min
./Allrun nonConformal    # the same with non-conformal periodic sides
./Allrun step 0.005      # afterwards: recording and rCFD run with another frame interval [s]
./Allrun reference       # afterwards: OpenFOAM solves the species on the same flow (7 min)
./Allrun clean
```

| Step | Files | Result |
|---|---|---|
| mesh | `makeLatticeMesh.py` | `constant/polyMesh` |
| steady flow | `system.flow`, `0/` | converged U, p, phi (`log.foamRun.flow`) |
| recording | `system.record`, `rCFD/rCFD_user_prep.h` | `data/`: two frames of the frozen flow |
| rCFD run | `rCFD/rCFD_user_run.h` | `rCFD_lattice_monitor.out`, `rCFD_lattice_profile.out`, `rCFD_results/` |
| reference | `system.reference`, `reference/` | `referenceRun/reference_monitor.out` |

The recording does not solve the flow again: `solver functions; subSolver incompressibleFluid;`
in `system.record/controlDict` reads the converged fields and runs only the function object.
Two frames are enough for a steady flow; the run replays them.

## The run header

One rCFD step:

1. `user_access_data_before_shift`: inlet cells set to A = 1.
2. Shift: every cell takes the content of its source cells of one frame interval earlier.
3. `user_access_data_after_swap`, function `diffuse_and_react`: in sub-steps, diffusion across
   all faces (`D |A| / d` per face, `rcfd_face_distance`, `rcfd_update_corona`) and the
   reactions in every cell. The wall reaction is in the cells that touch `reactingWall`
   (`rcfd_patch_faces`). Diffusion and wall reaction must be in the same sub-steps: the wall
   reaction is limited by the diffusion to the wall.
4. `user_post`: monitors and snapshots.

Settings at the top of `rCFD/rCFD_user_run.h`: `_global_time_step` (the frame interval),
`_diffusivity`, the rate constants, `_max_diffusion_number` (0.25: the sub-step follows from the
smallest cells), `_reaction_substep`, `_flow_rate`.

## Numbers

Mean mass fractions in the cells at the outlet after 4 s, 4 processors:

| Frame interval | A | C | D | E | Time of the run |
|---|---|---|---|---|---|
| reference (OpenFOAM, 2 ms steps) | 0.0052 | 0.1752 | 0.3404 | 0.4765 | 392 s |
| 100 ms | 0.0086 | 0.1653 | 0.3571 | 0.4665 | 27 s |
| 50 ms | 0.0103 | 0.1729 | 0.3463 | 0.4683 | 33 s |
| 20 ms | 0.0077 | 0.1727 | 0.3469 | 0.4702 | 43 s |
| 10 ms | 0.0041 | 0.1652 | 0.3546 | 0.4725 | 62 s |
| 5 ms | 0.0029 | 0.1583 | 0.3657 | 0.4678 | 88 s |

The products are within 0.03 of the reference for all frame intervals; 20 to 50 ms is the best
range here (the liquid moves 1 to 4 cells per frame in the lattice). Shorter frames do not
improve the result: the shifts carry a little less than the displacement when the origin of a
cell lies close to the cell itself (interpolation between the source cell and its face
neighbours), and that error adds up over more steps.

Cyclic sides and non-conformal sides give the same result (outlet values equal to four digits).

## A case of your own

1. Converged steady flow in OpenFOAM-13 (`incompressibleFluid`). A mesh with `cyclicAMI` sides
   from another OpenFOAM version: set the side patches to `type patch` in
   `constant/polyMesh/boundary` and couple them with `createNonConformalCouples`
   (`system.flow/createNonConformalCouplesDict`; `separation` is the position of the first patch
   minus the position of the second).
2. `system.record` and `rCFD/` from here. Patch names in `rCFD/rCFD_user_run.h`
   (`reactingWall`, `inlet`), flow rate, density, rate constants, `_z_min`/`_cell_size` of the
   profile output.
3. Frame interval so that the liquid moves one to a few cells per frame in the region of
   interest.
4. Cost to expect. Recording: the library searches the source cells by brute force inside each
   processor, minutes per frame for 250 000 cells per processor; two frames are needed. Memory:
   the coupling holds tables over all cells on every processor for a short time, about 75 bytes
   per cell of the whole mesh (1.2 GB per process for 16 million cells). Run: the diffusion
   sub-step scales with the square of the smallest cell size.
