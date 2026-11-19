GF25 fluidised bed, Euler-Euler (multiphaseEuler), OpenFOAM-13, with rCFD recording
Version 2 (new mesh, 15 kg bed, CFD-DEM parameters). Case files with geometry, without mesh.

Tested in this form on 4 processes for the first 20 ms: meshing, start-up, two frames of 5 ms,
closing of the database, three rcfd_run steps. Bed mass from the monitor: 15.01 kg.
The fluidised state itself has not been simulated yet.

Mesh (system/blockMeshDict, system/snappyHexMeshDict): 177 333 cells
  bed (y 0.04 to 0.43) 10 mm, nozzles and spray cones 5 mm, freeboard 20 mm (10 mm at the casing)
  largest non-orthogonality 48 degrees. Geometry: constant/geometry (six .stl and .eMesh files).

Taken from the CFD-DEM case:
  bed mass 15 kg                    system/setFieldsDict (alpha 0.5266 in 0.11 m)
  diameter 1.8 mm, density 1040     constant/phaseProperties, physicalProperties.particles
  drag Beetstra                     constant/momentumTransfer (BeetstraLimited, see dragModel/)
  restitution 0.90 / wall 0.75      constant/momentumTransport.particles, 0/Theta.particles
Still assumed:
  gas at 293.15 K, 1e5 Pa (rho 1.19 kg/m3, mu 1.84e-5 Pa s)     0/T.air, 0/T.particles
  wall specularity 0.01             0/U.particles, 0/Theta.particles
  nozzle velocity scaled to 20.67 m/s to keep the flow rate     BCDict

Use (run_gf25.sbatch does steps 2 to 5 and skips those already done):
  1. copy this folder to the work directory of the cluster
  2. of13;  (cd dragModel && wmake libso)      once: libGF25Drag.so
  3. blockMesh;  snappyHexMesh                  about one minute, serial
  4. rcfd_compile
  5. 32 processes: decomposePar; mpirun setFields -parallel; mpirun foamRun -parallel

Recording: 2000 frames of 5 ms from t = 10 s (rCFD/rCFD_user_prep.h, system/functions), endTime 20 s.
Averages (fieldAverage) start at 5 s.
Database: about 12 MB per frame, 25 GB in total.
Monitors in postProcessing/: solids volume (x 1040 = kg), pressure at the inlet, largest alpha.particles.
rCFD/rCFD_user_run.h is the package template (density 1040, 2000 frames, 5 ms): put your own model there.
run_gf25.sbatch: job script of the TUHH cluster (folder in its --chdir line).

Batch granulation (layering) on the same database: rCFD/rCFD_user_run_granulation.h
  Three moment fields of the particle diameter instead of size bins; the bed mass follows from the
  fields. Copy it to rCFD/rCFD_user_run.h (keep the mixing header under another name), rcfd_compile run.
  Settings at its top: seeds, densities of core and layer, spray rate, active nozzles, _growth_acceleration.
  Output: rCFD_granulation_monitor.out, rCFD_granulation_chambers.out; picture: py/rcfd_granulation_plot.py.
  Explanation with references: https://claude.ai/artifact/HfDBYMrcsFfC54utS7C8UH
  Tested on a small fluidised-bed database and for 40 steps on a two-frame GF25 database; a full batch
  on the recorded GF25 flow has not been run.

Batch in stages on several databases (bed mass grows beyond what one recorded flow stands for): run_gf25_batch.sbatch
  Databases of the same mesh at 15, 20, 25, ... kg (copy constant/polyMesh, do not mesh again).
  Batch folder with two files: rCFD_user_run.h (copy of rCFD/rCFD_user_run_granulation.h with your settings)
  and stages.txt (per line: case folder of a database, bed mass at which its stage ends in kg; the last
  value is the bed mass at the end of the batch). The job makes stage_1, stage_2, ... there, each with
  links to data/ and constant/ of its database, and writes rCFD_granulation.png for the whole batch.
  Needs package r6 (rcfd_run -init). Every stage must end within one job; with _growth_acceleration 10
  a batch of two hours needs about 8 hours on 32 processes.
