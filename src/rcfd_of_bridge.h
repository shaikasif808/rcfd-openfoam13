/* rcfd_of_bridge.h - rCFD preparation mode driven from a running solver (TUHH 11/2026)
 *
 * The solver (OpenFOAM function object rCFDPrep) hands over the mesh of its own partition
 * once and then one frame per frame interval. The dictionary settings come from the user
 * header rCFD_user_prep.h of the case, compiled into a shared library, exactly as with
 * the Fluent UDF. The bridge drives the unchanged rCFD library in preparation mode and
 * writes the database
 * (./data/c2c, vof, rec, rcfd_db.meta.toml) and, if asked, ./data/topo for the run mode.
 * All paths are relative to the current directory (the case directory).
 *
 * Floating point traps of the calling program (OpenFOAM: FOAM_SIGFPE) are switched off
 * inside the bridge and restored on return.
 *
 * MPI: the calling program's MPI_COMM_WORLD. If MPI is not initialised (serial solver
 * run), the bridge initialises it and finalises it in rcfd_of_prep_close().
 * Errors print a message and abort the program.
 */
#ifndef RCFD_OF_BRIDGE_H
#define RCFD_OF_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

/* the shared library is built with hidden visibility; only these functions are exported */
#if defined(__GNUC__)
#define RCFD_OF_API __attribute__((visibility("default")))
#else
#define RCFD_OF_API
#endif

typedef struct {

    int          n_int;             /* cells of this rank                                        */
    int          n_ext;             /* corona cells: cells of other ranks sharing a face         */
    int          n_fint;            /* faces between two cells of this rank                      */
    int          n_fext;            /* partition faces                                           */
    int          n_cells_global;    /* cells of the whole mesh                                   */

    /* cells: this rank's cells first, then the corona cells */
    const double *x;                /* [3 * n_cells] centres; a corona copy = the owner's values */
    const double *vol;              /* [n_cells]                                                 */
    const int    *host;             /* [n_cells] owning rank                                     */
    const int    *gid;              /* [n_cells] cell index in the undecomposed mesh             */

    /* faces: faces between own cells first, then partition faces (c0 own, c1 corona) */
    const int    *c0, *c1;          /* [n_faces] local cell indices                              */
    const double *area;             /* [3 * n_faces] area vector from c0 to c1                   */

    /* image cells behind periodic patches: the last n_img of the n_int cells of this rank.
     * An image is a copy of the cell with the same gid on the rank img_host, placed where
     * that cell would be if the mesh continued. The library treats images as cells of this
     * rank (neighbours for the interpolation of shifts); the bridge takes them out of the
     * database again. The faces to image cells are the last n_fimg of the n_fint faces:
     * one per pair of a cell and an image (across a non-conformal patch a cell can have
     * several images and an image several cells). */
    int          n_img;
    const int    *img_host;         /* [n_img] rank that holds the cell the image copies         */
    int          n_fimg;            /* faces between own cells and image cells                   */

} rcfd_of_mesh_t;

#define RCFD_OF_MAX_PHASES 8

typedef struct {

    const char *user_library;       /* shared library built from rCFD_user_prep.h (rcfd_compile)    */
    int     n_phases;               /* phases the solver hands over; must equal number_of_phases     */
    const char *phase_names;        /* "air particles": written to data/rcfd_foam_prep.txt           */
    int     norm_phase;             /* phase whose volume fraction is the recurrence norm [last]     */
    int     bbox_nodes;             /* [1]   padded bounding box of each rank                        */
    int     write_topo;             /* [1]   write data/topo/{topo_<rank>,gid_<rank>.bin,meta.txt}   */

} rcfd_of_prep_settings_t;

/* what the user header set; read back from the library after its initialisation */
typedef struct {

    int     n_frames;               /* Solver_Dict number_of_frames                                  */
    double  frame_interval;         /* Solver_Dict global_time_step                                  */
    int     n_phases;
    double  vof_cutoff;

} rcfd_of_prep_info_t;

RCFD_OF_API void rcfd_of_prep_default_settings(rcfd_of_prep_settings_t *s);

/* mesh arrays are copied; the user library is loaded and its user_register_functions()
 * called, as the Fluent UDF does with the user header; returns 0 */
RCFD_OF_API int  rcfd_of_prep_open(const rcfd_of_mesh_t *mesh, const rcfd_of_prep_settings_t *s, rcfd_of_prep_info_t *info);

/* one frame; arrays over the n_int + n_ext cells of this rank:
 *   vof[n_cells * n_phases]              (cell-major: vof[cell * n_phases + phase])
 *   origin[n_phases * 3 * n_cells]       (phase-major: origin[(phase * n_cells + cell) * 3 + d])
 * origin = position where the material that is in the cell now was one frame interval
 * earlier. Returns the number of frames written so far. */
RCFD_OF_API int  rcfd_of_prep_frame(const double *vof, const double *origin);
/* as rcfd_of_prep_frame; displ (same layout as origin, or NULL) is the vector from the cell
 * to its origin before the origin was moved through a periodic patch. It is used for the
 * displacement reported per frame only. */
RCFD_OF_API int  rcfd_of_prep_frame_displ(const double *vof, const double *origin, const double *displ);

RCFD_OF_API int  rcfd_of_prep_frames_written(void);

/* recurrence path, meta file, checks; returns 0 if the database is complete */
RCFD_OF_API int  rcfd_of_prep_close(void);

/* called from the user library (udf.h of this package) */
RCFD_OF_API int  rcfd_of_register_user_function(void *handle, const char *name, void (*f)(void));
RCFD_OF_API void rcfd_of_set_cell_dict(void *handle, const short i_layer);
RCFD_OF_API void rcfd_of_set_face_dict(void *handle, const short i_layer);

#ifdef __cplusplus
}
#endif

#endif
