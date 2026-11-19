/* udf.h - stand-in for the Fluent header when an rCFD user header is compiled for
 * OpenFOAM (rcfd_compile). TUHH 11/2026
 *
 * A user header written for Fluent keeps its first lines
 *
 *     #include "udf.h"
 *     #include "rCFD_library.h"
 *
 * and its functions (user_set_solver_dict, user_set_phase_dict, user_set_data_dict,
 * user_init_data, user_post, ..., user_register_functions). What Fluent offered and a
 * user header typically uses is defined here. Code that walks Fluent threads
 * (thread_loop_c, C_UDMI, ...) has no counterpart: put it between
 *
 *     #ifndef RCFD_OPENFOAM ... #endif
 *
 * and use rcfd_user_field() below to show per-cell values.
 */
#ifndef RCFD_FOAM_UDF_H
#define RCFD_FOAM_UDF_H

#define RCFD_OPENFOAM 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/* every process is a compute node */
#define RP_NODE  1
#define RP_HOST  0
#define PARALLEL 1

typedef double real;

extern int myid;                    /* MPI rank; set before user_register_functions is called */

#define node_zero         0
#define I_AM_NODE_ZERO_P  (myid == 0)

#define Message(...)   do{ printf(__VA_ARGS__); fflush(stdout); }while(0)
#define Message0(...)  do{ if(myid == 0){ printf(__VA_ARGS__); fflush(stdout); } }while(0)

/* global reductions over all ranks */
static inline double PRF_GRSUM1(double x)  { double g = x; MPI_Allreduce(&x, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return g; }
static inline double PRF_GRHIGH1(double x) { double g = x; MPI_Allreduce(&x, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return g; }
static inline double PRF_GRLOW1(double x)  { double g = x; MPI_Allreduce(&x, &g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); return g; }
static inline int    PRF_GISUM1(int x)     { int g = x;    MPI_Allreduce(&x, &g, 1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD); return g; }
static inline int    PRF_GIHIGH1(int x)    { int g = x;    MPI_Allreduce(&x, &g, 1, MPI_INT,    MPI_MAX, MPI_COMM_WORLD); return g; }
static inline int    PRF_GILOW1(int x)     { int g = x;    MPI_Allreduce(&x, &g, 1, MPI_INT,    MPI_MIN, MPI_COMM_WORLD); return g; }

/* registration of the user functions, used inside user_register_functions(void *handle) */
int rcfd_foam_register(void *handle, const char *name, void (*f)(void));

#define REGISTER_RCFD_UDF(funcName)  rcfd_foam_register(handle, #funcName, (void(*)(void))funcName)

/* Cell and face counts come from the OpenFOAM mesh (preparation) or from data/topo (run).
 * The Fluent calls are kept so that user_set_cell_dict / user_set_face_dict stay as they are. */
void fluent_set_cell_dict(void *handle, const short i_layer);
void fluent_set_face_dict(void *handle, const short i_layer);

/* Counterpart of a user-defined memory (C_UDMI) in run mode: a per-cell array with a name.
 * The first call creates it (length: internal cells of this rank, zero), later calls return
 * the same array. The run program writes it with every snapshot; snap_to_foam.py turns it
 * into the field rCFD_<name>. In preparation mode it returns NULL.
 *
 *     double *vf = rcfd_user_field(handle, "vf_diam_0");
 *     if(vf) vf[i_cell] = data_solid[i_cell][0] * cell_vof[i_cell][solid];
 */
double *rcfd_user_field(void *handle, const char *name);

/* Counterpart of a face loop over a boundary thread (begin_f_loop, F_C0) in run mode: the faces
 * of a boundary patch of the OpenFOAM mesh that lie on this processor.
 * Returns the number of faces here (0 if this processor has none) and -1 if the database has
 * no patch of that name (the log lists the names). cells[i] is the cell next to face i (an index
 * into the cell arrays), areas[i] the face area [m2], centres[3*i .. 3*i+2] the face centre [m].
 * Any of the three pointers may be NULL. In preparation mode it returns -1.
 *
 *     const int *c; int i, n = rcfd_patch_faces(handle, "inlet", &c, NULL, NULL);
 *     for(i = 0; i < n; i++) data_gas[c[i]][tracer] = 1.0;        // fixed value at an inlet
 */
int rcfd_patch_faces(void *handle, const char *patch, const int **cells, const double **areas, const double **centres);

/* For an exchange between neighbour cells written in the user header (a diffusion with a
 * diffusivity in m2/s, coupled with reactions in sub-steps), in run mode:
 *
 * rcfd_face_distance: distance[i] of the two cell centres of face i [m], for all faces of
 * rcfd_get_face_c0 / rcfd_get_face_c1. Across a periodic patch it is the distance through the
 * patch. Returns the number of faces, -1 in preparation mode.
 *
 * rcfd_update_corona: copies the values of all data of a phase from the cells to their copies
 * on the other processors (the cells c1 >= rcfd_get_int_cell_count of the faces at a processor
 * boundary). Call it before every loop over the faces. Returns 0.
 *
 *     const double *d; int n = rcfd_face_distance(handle, &d);
 *     rcfd_update_corona(handle, fluid);
 *     for(i = 0; i < n; i++){ exchange = D * area_i / d[i] * dt * (data[c1[i]][k] - data[c0[i]][k]); ... }
 */
int rcfd_face_distance(void *handle, const double **distance);
int rcfd_update_corona(void *handle, int i_phase);

/* rcfd_update_corona_field: the same for one array of the user header with one value per cell,
 * length rcfd_get_cell_count (cells of this processor first, then the copies of cells of the
 * other processors): the copies are filled. For iterative solvers written in the header (an
 * implicit diffusion step needs the search direction in the neighbour cells). Returns 0. */
int rcfd_update_corona_field(void *handle, double *field);

/* Runs in stages on several databases of the same mesh (run mode).
 *   mpirun -np N rCFD/rcfd_run ... -init DIR   starts with the data of every cell as they were at the end
 *   of the run in DIR (its out/state_final.bin, written at the end of every run).
 * rcfd_init_folder(): DIR, or NULL in a run without -init (available from user_init_data on).
 * rcfd_init_vof(handle, i_phase): volume fraction of the phase in every cell of this processor at the
 *   end of the run in DIR (the other database), NULL without -init; available from the first step on.
 *   A cell that was empty there holds data without meaning: the header has to set them.
 * rcfd_init_time(): time at the end of the run in DIR [s], 0 without -init.
 * rcfd_request_stop(handle): the run ends after the current step (call it on every processor or on one). */
const char   *rcfd_init_folder(void);
const double *rcfd_init_vof(void *handle, int i_phase);
double rcfd_init_time(void);
void   rcfd_request_stop(void *handle);

#ifdef __cplusplus
}
#endif

#endif
