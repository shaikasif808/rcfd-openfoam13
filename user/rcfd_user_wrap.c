/* rcfd_user_wrap.c - compiles the user header of a case (TUHH 11/2026)
 *
 * rcfd_compile builds this file twice:
 *   -DRCFD_USER_PREP -DRCFD_USER_HEADER='"rCFD_user_prep.h"'  -> rCFD/librCFD_user_prep.so,
 *        loaded by the function object rCFDPrep during the OpenFOAM run
 *   -DRCFD_USER_RUN  -DRCFD_USER_HEADER='"rCFD_user_run.h"'   -> linked into rCFD/rcfd_run
 * In both cases user_register_functions(handle) of the header is called before the
 * library is initialised, as the Fluent UDF does.
 */
#include "udf.h"
#include "rCFD_library.h"

int myid = 0;

#include RCFD_USER_HEADER

#if defined(RCFD_USER_PREP)

/* provided by librcfdOF.so (rcfd_of_bridge.c) */
int  rcfd_of_register_user_function(void *handle, const char *name, void (*f)(void));
void rcfd_of_set_cell_dict(void *handle, const short i_layer);
void rcfd_of_set_face_dict(void *handle, const short i_layer);

int    rcfd_foam_register(void *handle, const char *name, void (*f)(void)) { return rcfd_of_register_user_function(handle, name, f); }
void   fluent_set_cell_dict(void *handle, const short i_layer)             { rcfd_of_set_cell_dict(handle, i_layer); }
void   fluent_set_face_dict(void *handle, const short i_layer)             { rcfd_of_set_face_dict(handle, i_layer); }
double *rcfd_user_field(void *handle, const char *name)                    { (void)handle; (void)name; return NULL; }
int    rcfd_patch_faces(void *handle, const char *patch, const int **cells, const double **areas, const double **centres)
{
    (void)handle; (void)patch;

    if(cells) *cells = NULL;
    if(areas) *areas = NULL;
    if(centres) *centres = NULL;

    return -1;      /* run mode only */
}

#else

/* provided by the run program (rcfd_run_main.c) */
int    rcfd_run_register_user_function(void *handle, const char *name, void (*f)(void));
double *rcfd_run_user_field(void *handle, const char *name);
int    rcfd_run_patch_faces(void *handle, const char *patch, const int **cells, const double **areas, const double **centres);

int    rcfd_foam_register(void *handle, const char *name, void (*f)(void)) { return rcfd_run_register_user_function(handle, name, f); }
void   fluent_set_cell_dict(void *handle, const short i_layer)             { (void)handle; (void)i_layer; }   /* counts are read from data/topo */
void   fluent_set_face_dict(void *handle, const short i_layer)             { (void)handle; (void)i_layer; }
double *rcfd_user_field(void *handle, const char *name)                    { return rcfd_run_user_field(handle, name); }
int    rcfd_patch_faces(void *handle, const char *patch, const int **cells, const double **areas, const double **centres)
{
    return rcfd_run_patch_faces(handle, patch, cells, areas, centres);
}

#endif

/* cksum and size of the header at build time (set by rcfd_compile); the function object and
 * the run program stop if the header in the case was changed after the build */
#ifndef RCFD_USER_HEADER_CRC
#define RCFD_USER_HEADER_CRC    0UL
#define RCFD_USER_HEADER_BYTES  -1L
#endif

unsigned long rcfd_user_header_crc(void)   { return RCFD_USER_HEADER_CRC; }
long          rcfd_user_header_bytes(void) { return RCFD_USER_HEADER_BYTES; }

/* entry point: the same name in the preparation library and in the run program */
void rcfd_user_register_functions(void *handle)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    user_register_functions(handle);
}
