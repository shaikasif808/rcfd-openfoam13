fischerK_of13: the 16 million cell lattice case (fischerK_a1_0_b1_0_c1_0) in OpenFOAM-13 for rCFD

Job 1 (this package): mesh and steady flow.           sbatch job1_flow.sbatch
Job 2 (follows):      rCFD recording on the frozen flow and rCFD run with the reactions.

What job 1 does
  1. Copies points, faces, owner, neighbour and boundary from
     ../fischerK_a1_0_b1_0_c1_0/constant/polyMesh (3.6 GB). The v2512 case is only read.
     Zones and sets are not copied.
  2. boundary: the four cyclicAMI patches (left, right, top, bottom) become ordinary patches;
     createNonConformalCouples couples left with right and top with bottom
     (system/createNonConformalCouplesDict).
  3. decomposePar (64, scotch; hierarchical if the build has no scotch), checkMesh.
  4. foamRun: steady laminar flow, incompressibleFluid, 1 mL/s at the inlet,
     nu = 0.6527e-3 / 992.25 m2/s. Stops when the residuals of p and U are below 1e-5 and 1e-6
     (system/fvSolution) or after 3000 iterations.

Submitting again skips what is done and continues the flow from the last written iteration.
To start from scratch: remove constant/polyMesh, constant/fvMesh, processor*, 0, postProcessing
and log.* in this folder.

Disk space: about 10 GB in this folder.
Run time: not measured. Estimate from the small test case: copy, coupling and decomposition
30 to 60 min (serial), flow 1 to 3 s per iteration.

Tested with the small lattice mesh (209 104 cells) written in the file format of the v2512
mesh, 4 processes: conversion, coupling (coverage 1), checkMesh, flow converged in 241
iterations, continuation after a second start. Not tested: the 16 million cell mesh itself.

What to send back: the file log.slurm.<job number>.
