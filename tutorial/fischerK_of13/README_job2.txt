fischerK_of13, job 2: rCFD database on the frozen flow and rCFD run.      sbatch job2_rcfd.sbatch

Needs: job 1 finished (processor*/<last iteration> with phi), package r5 of rcfd_foam2 built.

Files of this package
  system.record/        recording: solver functions on the converged flow, function object
                        rCFDPrep, frames of 0.01 s, time step 0.5 ms
  rCFD/rCFD_user_prep.h frame interval 0.01 s, two frames
  rCFD/rCFD_user_run.h  A at the inlet, wall reaction at reactingWall, B -> C, C -> D, 2 B -> E,
                        diffusion D = nu; implicit diffusion step with the wall reaction,
                        4 sub-steps per rCFD step (_diffusion_in_header 2, _implicit_substeps)
  job2_rcfd.sbatch      recording, rCFD run (500 steps = 5 s), snapshots for ParaView

What the job changes in the folder
  processor*/0          replaced by a copy of the converged flow (the flow result stays)
  system                renamed to system.flow; system.record is copied to system
  data/                 the database (new)
  out/, rCFD_*.out, rCFD_results/, log.*          results (earlier ones are moved to rcfd_run_<date>/)

Output
  rCFD_lattice_monitor.out   per step: time, mean A B C D E in the cells at the outlet,
                             content of the domain [kg]
  rCFD_lattice_profile.out   volume-weighted mean A B C D E in 70 layers of 0.5 mm along z
                             (z = -10 ... 25 mm): compare with the plane averages of the
                             reacting case (speciesAverages, totalFlux.outlet at z = 10 mm)
  log.rcfd_run               "[user] step N: iterations of the diffusion step ..."

Tested with the small lattice mesh (209 104 cells, non-conformal sides, 4 processes): outlet
after 4 s A 0.0040, C 0.1654, D 0.3570, E 0.4698; OpenFOAM on the same flow 0.0052, 0.1752,
0.3404, 0.4765. Not tested: the 16 million cell mesh. Unknown there: the time of the recording
(estimate 10 to 20 min) and the number of iterations of the diffusion step (it decides the
time of the run). A short test first:   sbatch --export=ALL,STEPS=20,SNAP=0 job2_rcfd.sbatch
