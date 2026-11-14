/* rcfd_of_bridge.c - rCFD preparation mode driven from a running solver (TUHH 11/2026)
 *
 * See rcfd_of_bridge.h. Mesh and frames come through memory from the solver. The
 * dictionary settings come from the user header of the case (rCFD_user_prep.h), which
 * rcfd_compile builds into a shared library; it is loaded here and its
 * user_register_functions() is called, as the Fluent UDF does.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include <fenv.h>
#include <dlfcn.h>

#include "rcfd_of_common.h"
#include "rcfd_of_bridge.h"

static void     *g_h = NULL;
static int      g_open = 0, g_mpi_started_here = 0;
static int      g_n_frames = 0, g_frames_written = 0, g_n_phases = 0, g_norm_phase = 0;
static int      g_norm_count = 0, g_norm_coarse = 1;
static double   g_cutoff = 0.0, g_frame_interval = 0.0;
static double   *g_displ = NULL, *g_norms = NULL, **g_vof = NULL, *g_vof_mem = NULL;
static double   g_t_shift = 0.0, g_t_interp = 0.0, g_t_write = 0.0;
static long long g_total_shifts[RCFD_OF_MAX_PHASES], g_total_incoming[RCFD_OF_MAX_PHASES];
static long long g_c2c_offset[RCFD_OF_MAX_PHASES];
static rcfd_of_prep_settings_t g_set;
static void     *g_user_lib = NULL;
static int      g_interpolation = -1;      /* -1 not known yet, 0 off, 1 on */

/* ------------------------------------------------------------------------- */
/* user functions of the case                                                 */
/*                                                                            */
/* The user header registers its functions through REGISTER_RCFD_UDF, which   */
/* ends up in rcfd_of_register_user_function(). Everything goes straight to   */
/* the library except the solver dict: after the user's function the bridge   */
/* sets the two properties the OpenFOAM coupling depends on.                  */
/* ------------------------------------------------------------------------- */

static void (*g_user_solver_dict)(void *) = NULL;

static void bridge_set_solver_dict(void *h)
{
    if(g_user_solver_dict) g_user_solver_dict(h);

    /* the function object hands over position labels, one frame per call */
    rcfd_set_dict_property_i(h, "Solver_Dict", -1, -1, -1, "shifts_from_displacement_fields_on", 1);
    rcfd_set_dict_property_i(h, "Solver_Dict", -1, -1, -1, "time_steps_per_monitoring_interval", 1);
}

int rcfd_of_register_user_function(void *handle, const char *name, void (*f)(void))
{
    if(strcmp(name, "user_set_solver_dict") == 0){

        g_user_solver_dict = (void(*)(void*))f;

        return 0;       /* bridge_set_solver_dict is registered and calls it */
    }

    if(strcmp(name, "user_init_cells") == 0 || strcmp(name, "user_init_faces") == 0 ||
       strcmp(name, "user_init_cell_hosts") == 0 || strcmp(name, "user_init_nodes") == 0 ||
       strcmp(name, "user_set_node_dict") == 0){

        of_msg("note: %s of the user header is not used; cells, faces and hosts come from the OpenFOAM mesh", name);

        return 0;
    }

    if(rcfd_register_user_function(handle, name, f) != 0){

        of_msg("WARNING: the library does not know the user function '%s' (ignored)", name);

        return -1;
    }

    return 0;
}

/* fluent_set_cell_dict / fluent_set_face_dict of the Fluent headers map to these */
void rcfd_of_set_cell_dict(void *handle, const short i_layer) { of_user_set_cell_dict(handle, i_layer); }
void rcfd_of_set_face_dict(void *handle, const short i_layer) { of_user_set_face_dict(handle, i_layer); }

static void bridge_load_user_library(void *h, const char *path)
{
    void (*reg)(void *);
    char full[4096];

    if(path == NULL || path[0] == 0) of_abort("no user library given (function object entry userLibrary)");

    if(path[0] == '/') snprintf(full, sizeof(full), "%s", path);
    else               snprintf(full, sizeof(full), "./%s", path);

    if(!of_file_exists(full))
        of_abort("user library '%s' not found. Put rCFD_user_prep.h into the folder rCFD of the case and run rcfd_compile there", full);

    g_user_lib = dlopen(full, RTLD_NOW | RTLD_LOCAL);

    if(g_user_lib == NULL) of_abort("cannot load the user library '%s': %s", full, dlerror());

    reg = (void(*)(void*))dlsym(g_user_lib, "rcfd_user_register_functions");

    if(reg == NULL) of_abort("'%s' has no rcfd_user_register_functions: build it with rcfd_compile", full);

    /* rCFD_user_prep.h next to the library: stop if it was changed after the build */
    {
        unsigned long (*crc)(void) = (unsigned long(*)(void))dlsym(g_user_lib, "rcfd_user_header_crc");
        long (*bytes)(void)        = (long(*)(void))dlsym(g_user_lib, "rcfd_user_header_bytes");
        char header[4200], *slash = strrchr(full, '/');

        if(crc != NULL && bytes != NULL && slash != NULL){

            snprintf(header, sizeof(header), "%.*s/rCFD_user_prep.h", (int)(slash - full), full);

            of_check_user_header(header, full, crc(), bytes(), "rcfd_compile prep");
        }
    }

    reg(h);     /* -> user_register_functions(handle) of rCFD_user_prep.h */
}

/* ------------------------------------------------------------------------- */
/* number of shifts the library wrote for the last frame                      */
/*                                                                            */
/* binary c2c file, per frame: int n_interior, int n_incoming, then           */
/* n_interior x (int c0, int c1, uchar w) and n_incoming x (int c0, node0,    */
/* int c1, uchar w), node0 one byte up to 256 ranks, else int                 */
/* ------------------------------------------------------------------------- */

static void prep_count_written_shifts(int i_phase, long long *n_interior, long long *n_incoming)
{
    char  name[256];
    FILE  *f;
    int   head[2] = {0, 0};
    long long rec_in = 9, rec_out = (of_size > 256) ? 13 : 10, size;

    *n_interior = 0; *n_incoming = 0;

    snprintf(name, sizeof(name), "data/c2c/c2c_0_%d_%d", i_phase, of_rank);

    size = of_file_size(name);

    if(size < g_c2c_offset[i_phase] + 8) of_abort("%s: the library did not write the shifts of this frame", name);

    f = fopen(name, "rb");

    if(f == NULL) of_abort("cannot open '%s'", name);

    if(fseek(f, (long)g_c2c_offset[i_phase], SEEK_SET) != 0 || fread(head, sizeof(int), 2, f) != 2) of_abort("%s: cannot read the frame header", name);

    fclose(f);

    *n_interior = head[0];
    *n_incoming = head[1];

    g_c2c_offset[i_phase] += 8 + rec_in * head[0] + rec_out * head[1];

    if(g_c2c_offset[i_phase] != size)
        of_abort("%s: size %lld after this frame, expected %lld from the frame header (%d interior, %d incoming shifts)",
                 name, size, g_c2c_offset[i_phase], head[0], head[1]);
}

/* remove database files of ranks >= of_size left over from a preparation on more ranks */
static void prep_remove_stale_files(void)
{
    static const char *pattern[5] = {"data/c2c/c2c_0_%d_%d", "data/vof/vof_0_%d_%d", "data/tmp/norm_0_%d", "data/topo/topo_%d", "data/topo/gid_%d.bin"};
    char name[256];
    int  r, p, k, removed = 0, found;

    if(of_rank != 0) return;

    for(r = of_size; r < of_size + 4096; r++){

        found = 0;

        for(k = 0; k < 5; k++){

            if(k >= 3 && !g_set.write_topo) continue;

            for(p = 0; p < RCFD_OF_MAX_PHASES; p++){

                if(k >= 2){ if(p > 0) break; snprintf(name, sizeof(name), pattern[k], r); }
                else       { snprintf(name, sizeof(name), pattern[k], p, r); }

                if(of_file_exists(name)){ found = 1; if(remove(name) == 0) removed++; }
            }
        }

        if(!found) break;
    }

    if(removed > 0) of_msg("removed %d database files of ranks >= %d (left over from an earlier preparation)", removed, of_size);
}

/* ------------------------------------------------------------------------- */
/* mesh from the caller                                                       */
/* ------------------------------------------------------------------------- */

static void bridge_copy_mesh(of_topo_t *T, const rcfd_of_mesh_t *m)
{
    int i, n_cells, n_faces;

    memset(T, 0, sizeof(*T));

    if(m->n_int < 1)   of_abort("this rank has no cells (a rank without cells is not supported by the library)");
    if(m->n_ext < 0 || m->n_fint < 0 || m->n_fext < 0) of_abort("negative count in the mesh handed to the bridge");
    if(m->n_fint + m->n_fext < 1) of_abort("this rank has no faces (a rank without faces is not supported by the library)");

    n_cells = m->n_int + m->n_ext;
    n_faces = m->n_fint + m->n_fext;

    T->n_cells = n_cells; T->n_int = m->n_int; T->n_ext = m->n_ext;
    T->n_faces = n_faces; T->n_fint = m->n_fint; T->n_fext = m->n_fext;
    T->n_nodes = 0;

    T->x    = (double*)of_malloc(3 * (size_t)n_cells * sizeof(double), "cell centres");
    T->vol  = (double*)of_malloc((size_t)n_cells * sizeof(double), "cell volumes");
    T->host = (int*)   of_malloc((size_t)n_cells * sizeof(int), "cell hosts");
    T->gid  = (int*)   of_malloc((size_t)n_cells * sizeof(int), "global cell ids");
    T->c0   = (int*)   of_malloc((size_t)n_faces * sizeof(int), "face c0");
    T->c1   = (int*)   of_malloc((size_t)n_faces * sizeof(int), "face c1");
    T->area = (double*)of_malloc(3 * (size_t)n_faces * sizeof(double), "face areas");

    memcpy(T->x,    m->x,    3 * (size_t)n_cells * sizeof(double));
    memcpy(T->vol,  m->vol,  (size_t)n_cells * sizeof(double));
    memcpy(T->host, m->host, (size_t)n_cells * sizeof(int));
    memcpy(T->gid,  m->gid,  (size_t)n_cells * sizeof(int));
    memcpy(T->c0,   m->c0,   (size_t)n_faces * sizeof(int));
    memcpy(T->c1,   m->c1,   (size_t)n_faces * sizeof(int));
    memcpy(T->area, m->area, 3 * (size_t)n_faces * sizeof(double));

    for(i = 0; i < n_cells; i++){

        if(T->host[i] < 0 || T->host[i] >= of_size) of_abort("cell %d has host %d, but there are %d ranks", i, T->host[i], of_size);

        if(i < T->n_int && T->host[i] != of_rank)  of_abort("cell %d is listed as a cell of this rank, but its host is %d", i, T->host[i]);
        if(i >= T->n_int && T->host[i] == of_rank) of_abort("corona cell %d has this rank as host", i);

        if(!(T->vol[i] > 0.0)) of_abort("cell %d has volume %g", i, T->vol[i]);

        if(T->gid[i] < 0 || T->gid[i] >= of_n_cells_global) of_abort("cell %d has global index %d, outside 0..%d", i, T->gid[i], of_n_cells_global - 1);
    }

    for(i = 0; i < n_faces; i++){

        if(T->c0[i] < 0 || T->c0[i] >= T->n_int) of_abort("face %d: c0 = %d is not a cell of this rank (0..%d)", i, T->c0[i], T->n_int - 1);

        if(i < T->n_fint){

            if(T->c1[i] < 0 || T->c1[i] >= T->n_int || T->c1[i] == T->c0[i]) of_abort("face %d: c1 = %d invalid for a face between own cells (c0 = %d)", i, T->c1[i], T->c0[i]);
        }
        else{

            if(T->c1[i] < T->n_int || T->c1[i] >= n_cells) of_abort("partition face %d: c1 = %d is not a corona cell (%d..%d)", i, T->c1[i], T->n_int, n_cells - 1);
        }
    }
}

/* a corona copy must carry exactly the centre and volume of the owner's cell, and every
 * partition face must exist on both ranks with opposite area vectors */
static void bridge_check_corona(const of_topo_t *T)
{
    double  *loc, *glob, sum_a[3] = {0.0, 0.0, 0.0}, g_sum_a[3], abs_a = 0.0, g_abs_a = 0.0;
    int     i, d, bad = 0, g_bad = 0, *used, unused = 0, g_unused = 0;
    long long n_pf = T->n_fext, g_n_pf = 0;
    size_t  n = (size_t)of_n_cells_global;

    loc  = (double*)calloc(4 * n, sizeof(double));
    glob = (double*)calloc(4 * n, sizeof(double));

    if(loc == NULL || glob == NULL) of_abort("out of memory (corona check)");

    for(i = 0; i < T->n_int; i++){

        for(d = 0; d < 3; d++) loc[4 * (size_t)T->gid[i] + d] = T->x[3 * i + d];

        loc[4 * (size_t)T->gid[i] + 3] = T->vol[i];
    }

    MPI_Allreduce(loc, glob, (int)(4 * n), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(i = T->n_int; i < T->n_cells; i++){

        const double *g = &glob[4 * (size_t)T->gid[i]];

        if(g[0] != T->x[3 * i] || g[1] != T->x[3 * i + 1] || g[2] != T->x[3 * i + 2] || g[3] != T->vol[i]) bad++;
    }

    free(loc); free(glob);

    /* every corona cell must be c1 of a partition face */
    used = (int*)calloc((size_t)T->n_cells, sizeof(int));

    if(used == NULL) of_abort("out of memory (corona check)");

    for(i = T->n_fint; i < T->n_faces; i++){

        used[T->c1[i]] = 1;

        for(d = 0; d < 3; d++){ sum_a[d] += T->area[3 * i + d]; abs_a += fabs(T->area[3 * i + d]); }
    }

    for(i = T->n_int; i < T->n_cells; i++){ if(!used[i]) unused++; }

    free(used);

    MPI_Allreduce(&bad, &g_bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&unused, &g_unused, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&n_pf, &g_n_pf, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(sum_a, g_sum_a, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&abs_a, &g_abs_a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(g_bad > 0)    of_abort("%d corona cells do not carry the centre and volume of their owner (all ranks)", g_bad);
    if(g_unused > 0) of_abort("%d corona cells are not connected by a partition face (all ranks)", g_unused);
    if(g_n_pf % 2)   of_abort("the partition faces of all ranks add up to an odd number (%lld); every one must be listed on two ranks", g_n_pf);

    if(g_abs_a > 0.0 && (fabs(g_sum_a[0]) + fabs(g_sum_a[1]) + fabs(g_sum_a[2])) > 1.0e-9 * g_abs_a)
        of_abort("the area vectors of the partition faces of all ranks do not cancel (sum %g %g %g): they must point from the own to the corona cell",
                 g_sum_a[0], g_sum_a[1], g_sum_a[2]);
}

static void bridge_mkdir(const char *name)
{
    if(mkdir(name, 0777) != 0 && errno != EEXIST) of_abort("cannot create directory '%s': %s", name, strerror(errno));
}

/* data/topo/topo_<rank> (the library's ASCII format, no node section), gid_<rank>.bin, meta.txt */
static void bridge_write_topo(const of_topo_t *T)
{
    char    name[256];
    FILE    *f;
    int     i;

    bridge_mkdir("data");
    bridge_mkdir("data/topo");

    snprintf(name, sizeof(name), "data/topo/topo_%d", of_rank);

    f = fopen(name, "w");

    if(f == NULL) of_abort("cannot write '%s': %s", name, strerror(errno));

    fprintf(f, "%d %d %d \n", T->n_cells, T->n_faces, 0);

    for(i = 0; i < T->n_cells; i++)
        fprintf(f, "%.17e %.17e %.17e %.17e %d \n", T->x[3 * i], T->x[3 * i + 1], T->x[3 * i + 2], T->vol[i], T->host[i]);

    for(i = 0; i < T->n_faces; i++)
        fprintf(f, "%d %d %.17e %.17e %.17e \n", T->c0[i], T->c1[i], T->area[3 * i], T->area[3 * i + 1], T->area[3 * i + 2]);

    if(fclose(f) != 0) of_abort("error while writing '%s'", name);

    snprintf(name, sizeof(name), "data/topo/gid_%d.bin", of_rank);

    f = fopen(name, "wb");

    if(f == NULL) of_abort("cannot write '%s': %s", name, strerror(errno));

    {
        int32_t *buf = (int32_t*)of_malloc((size_t)T->n_cells * sizeof(int32_t), "gid buffer");

        for(i = 0; i < T->n_cells; i++) buf[i] = (int32_t)T->gid[i];

        if(of_fwrite_le(buf, sizeof(int32_t), (size_t)T->n_cells, f) != (size_t)T->n_cells) of_abort("error while writing '%s'", name);

        free(buf);
    }

    if(fclose(f) != 0) of_abort("error while writing '%s'", name);

    /* cell_rank.bin: int32 rank of every cell of the whole mesh (for post-processing) */
    {
        int32_t *loc  = (int32_t*)calloc((size_t)of_n_cells_global, sizeof(int32_t));
        int32_t *glob = (int32_t*)calloc((size_t)of_n_cells_global, sizeof(int32_t));

        if(loc == NULL || glob == NULL) of_abort("out of memory (cell ranks)");

        for(i = 0; i < T->n_int; i++) loc[T->gid[i]] = (int32_t)of_rank;

        MPI_Allreduce(loc, glob, of_n_cells_global, MPI_INT32_T, MPI_SUM, MPI_COMM_WORLD);

        if(of_rank == 0){

            f = fopen("data/topo/cell_rank.bin", "wb");

            if(f == NULL) of_abort("cannot write 'data/topo/cell_rank.bin': %s", strerror(errno));

            if(of_fwrite_le(glob, sizeof(int32_t), (size_t)of_n_cells_global, f) != (size_t)of_n_cells_global) of_abort("error while writing 'data/topo/cell_rank.bin'");

            fclose(f);
        }

        free(loc); free(glob);
    }

    if(of_rank == 0){

        f = fopen("data/topo/meta.txt", "w");

        if(f == NULL) of_abort("cannot write 'data/topo/meta.txt': %s", strerror(errno));

        fprintf(f, "%d %d\n", of_n_cells_global, of_size);

        fclose(f);
    }
}

/* case.cfg for rcfd_of_run, only if there is none */

/* settings of the preparation for the run program: it compares them with its own header */
static void bridge_write_prep_info(void)
{
    FILE *f;

    if(of_rank != 0) return;

    f = fopen("data/rcfd_foam_prep.txt", "w");

    if(f == NULL){ of_msg("WARNING: cannot write data/rcfd_foam_prep.txt: %s", strerror(errno)); return; }

    fprintf(f, "frame_interval = %.17g\n", g_frame_interval);
    fprintf(f, "n_frames = %d\n", g_n_frames);
    fprintf(f, "n_phases = %d\n", g_n_phases);
    fprintf(f, "n_ranks = %d\n", of_size);
    fprintf(f, "phases = %s\n", g_set.phase_names ? g_set.phase_names : "");

    fclose(f);
}

/* ------------------------------------------------------------------------- */
/* public functions                                                           */
/* ------------------------------------------------------------------------- */

void rcfd_of_prep_default_settings(rcfd_of_prep_settings_t *s)
{
    s->user_library = "rCFD/librCFD_user_prep.so";
    s->n_phases     = 0;
    s->phase_names  = NULL;
    s->norm_phase   = -1;
    s->bbox_nodes   = 1;
    s->write_topo   = 1;
}

static int bridge_open(const rcfd_of_mesh_t *mesh, const rcfd_of_prep_settings_t *s, rcfd_of_prep_info_t *info)
{
    int     i, flag = 0, n_cells;
    void    *h = NULL;

    if(g_open){ fprintf(stderr, "[rcfd_of] ERROR: rcfd_of_prep_open called twice\n"); return -1; }

    MPI_Initialized(&flag);

    if(!flag){ MPI_Init(NULL, NULL); g_mpi_started_here = 1; }

    MPI_Comm_rank(MPI_COMM_WORLD, &of_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &of_size);

    (void)of_fread_le; (void)of_read_frames_meta; (void)of_read_cfg; (void)of_read_topo_meta; (void)of_read_topo; (void)of_startup;

    g_set = *s;

    if(s->n_phases < 1 || s->n_phases > RCFD_OF_MAX_PHASES) of_abort("%d phases handed over (1 to %d are possible)", s->n_phases, RCFD_OF_MAX_PHASES);

    of_cfg_defaults(&of_cfg);

    of_cfg.prep_bbox_nodes = s->bbox_nodes;

    g_frames_written   = 0;
    g_interpolation    = -1;
    g_user_solver_dict = NULL;
    of_n_cells_global  = mesh->n_cells_global;
    of_meta_ranks      = of_size;

    if(of_n_cells_global < 1) of_abort("number of cells of the whole mesh: %d", of_n_cells_global);

    bridge_copy_mesh(&of_topo, mesh);
    of_check_topo_global(&of_topo);
    bridge_check_corona(&of_topo);

    {
        long long loc[4], sum[4];

        loc[0] = of_topo.n_int; loc[1] = of_topo.n_ext; loc[2] = of_topo.n_fint; loc[3] = of_topo.n_fext;

        MPI_Allreduce(loc, sum, 4, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        of_msg("preparation: %d ranks, %d cells; all ranks: %lld corona cells, %lld internal faces, %lld partition faces (counted on both sides)",
               of_size, of_n_cells_global, sum[1], sum[2], sum[3]);
    }

    /* library in preparation mode: topology from the solver, everything else from the user header */
    rcfd_open(0, NULL, MPI_COMM_WORLD, &h);

    of_register_topology(h, 1);

    rcfd_register_user_function(h, "user_set_solver_dict", (void(*)(void))bridge_set_solver_dict);

    bridge_load_user_library(h, s->user_library);

    rcfd_init(h, RCFD_MODE_PREPARATION);

    if(!rcfd_is_initialized(h)) of_abort("rcfd_init failed");

    of_check_library_topology(h, &of_topo);

    g_n_frames       = rcfd_get_frame_count(h);
    g_n_phases       = rcfd_get_phase_count(h);
    g_frame_interval = rcfd_get_global_timestep(h);
    g_cutoff         = rcfd_get_vof_cutoff_limit(h);

    if(g_n_phases != s->n_phases)
        of_abort("rCFD_user_prep.h sets number_of_phases = %d, the function object hands over %d phases (%s)",
                 g_n_phases, s->n_phases, s->phase_names ? s->phase_names : "");

    if(g_n_frames < 2)             of_abort("rCFD_user_prep.h: number_of_frames = %d (at least 2)", g_n_frames);
    if(!(g_frame_interval > 0.0))  of_abort("rCFD_user_prep.h: global_time_step = %g s", g_frame_interval);

    if(!rcfd_is_shifts_from_displacement_fields_on(h)) of_abort("the library did not take shifts_from_displacement_fields_on");

    g_norm_phase = (s->norm_phase >= 0 && s->norm_phase < g_n_phases) ? s->norm_phase : g_n_phases - 1;

    of_msg("preparation (rCFD_user_prep.h): %d frames of %g s, %d phases (%s), vof cutoff %g",
           g_n_frames, g_frame_interval, g_n_phases, s->phase_names ? s->phase_names : "", g_cutoff);

    prep_remove_stale_files();

    if(s->write_topo) bridge_write_topo(&of_topo);

    bridge_mkdir("data");
    bridge_write_prep_info();

    n_cells       = of_topo.n_cells;
    g_norm_count  = rcfd_get_norm_count(h);
    g_norm_coarse = rcfd_get_norm_coarsening(h);

    if(g_norm_coarse < 1) g_norm_coarse = 1;

    if(of_cfg.prep_bbox_nodes) of_msg("bounding box of each rank = box of its cell centres padded by its largest cell size (rank 0: %.4g m)", of_bbox_pad);

    g_displ   = (double*) of_malloc(3 * (size_t)n_cells * sizeof(double), "displacement field");
    g_norms   = (double*) of_malloc((size_t)(g_norm_count > 0 ? g_norm_count : 1) * sizeof(double), "norms");
    g_vof_mem = (double*) of_malloc((size_t)g_n_phases * (size_t)n_cells * sizeof(double), "vof");
    g_vof     = (double**)of_malloc((size_t)n_cells * sizeof(double*), "vof rows");

    for(i = 0; i < n_cells; i++) g_vof[i] = &g_vof_mem[(size_t)g_n_phases * i];

    for(i = 0; i < RCFD_OF_MAX_PHASES; i++){ g_c2c_offset[i] = 0; g_total_shifts[i] = 0; g_total_incoming[i] = 0; }

    g_t_shift = 0.0; g_t_interp = 0.0; g_t_write = 0.0;

    if(info){

        info->n_frames       = g_n_frames;
        info->frame_interval = g_frame_interval;
        info->n_phases       = g_n_phases;
        info->vof_cutoff     = g_cutoff;
    }

    g_h    = h;
    g_open = 1;

    return 0;
}

int rcfd_of_prep_frames_written(void)
{
    return g_frames_written;
}

static int bridge_frame(const double *vof_in, const double *origin_in)
{
    void    *h = g_h;
    int     i, d, p, err, i_frame = g_frames_written, np = g_n_phases;
    int     n_cells = of_topo.n_cells, n_int = of_topo.n_int;
    long long loc[3 * RCFD_OF_MAX_PHASES], glob[3 * RCFD_OF_MAX_PHASES];
    double  t0, ts, ti, tw;
    double  dl[RCFD_OF_MAX_PHASES], dg[RCFD_OF_MAX_PHASES], ml[RCFD_OF_MAX_PHASES], mg[RCFD_OF_MAX_PHASES];

    if(!g_open) of_abort("rcfd_of_prep_frame called before rcfd_of_prep_open");

    if(i_frame >= g_n_frames) return g_frames_written;      /* database is full */

    t0 = MPI_Wtime();

    for(p = 0; p < 3 * RCFD_OF_MAX_PHASES; p++) loc[p] = 0;

    for(i = 0; i < n_cells; i++){

        for(p = 0; p < np; p++){

            double a = vof_in[(size_t)np * i + p];

            if(!(a >= -1.0e-6 && a <= 1.0 + 1.0e-6))
                of_abort("frame %d: volume fraction of phase %d = %g in global cell %d", i_frame, p, a, of_topo.gid[i]);

            if(a < 0.0) a = 0.0;
            if(a > 1.0) a = 1.0;

            g_vof[i][p] = a;

            if(i < n_int && a > g_cutoff) loc[p]++;
        }
    }

    if(rcfd_set_vof(h, 0, g_vof) != 0) of_abort("rcfd_set_vof failed (frame %d)", i_frame);

    for(p = 0; p < np; p++){

        /* the library calls these absolute positions cfd_displ. Cells of a phase at or
           below the vof cutoff get 0. */
        const double *origin = origin_in + (size_t)p * 3 * (size_t)n_cells;

        dl[p] = 0.0; ml[p] = 0.0;

        for(i = 0; i < n_cells; i++){

            if(g_vof[i][p] > g_cutoff){

                double s2 = 0.0;

                for(d = 0; d < 3; d++){

                    double v = origin[3 * (size_t)i + d];

                    if(!(v == v) || fabs(v) > 1.0e30)
                        of_abort("frame %d: origin of phase %d in global cell %d is %g", i_frame, p, of_topo.gid[i], v);

                    g_displ[3 * i + d] = v;

                    s2 += (v - of_topo.x[3 * i + d]) * (v - of_topo.x[3 * i + d]);
                }

                if(i < n_int){

                    dl[p] += s2;

                    if(s2 > ml[p]) ml[p] = s2;
                }
            }
            else{

                g_displ[3 * i] = 0.0; g_displ[3 * i + 1] = 0.0; g_displ[3 * i + 2] = 0.0;
            }
        }

        if(rcfd_set_displacement_fields(h, 0, p, g_displ, NULL) != 0)
            of_abort("rcfd_set_displacement_fields failed (frame %d, phase %d)", i_frame, p);
    }

    if(rcfd_create_shifts_from_displacement_fields(h) != 0)
        of_abort("rcfd_create_shifts_from_displacement_fields failed (frame %d)", i_frame);

    ts = MPI_Wtime() - t0;

    t0 = MPI_Wtime();

    /* the library has no getter for interpolation_on: it returns -1 at once if it is off */
    if(g_interpolation != 0){

        if(rcfd_interpolate_shifts(h) != 0){

            if(g_interpolation < 0){ g_interpolation = 0; of_msg("interpolated shifts are off (interpolation_on in rCFD_user_prep.h)"); }
            else                    of_abort("rcfd_interpolate_shifts failed (frame %d)", i_frame);
        }
        else g_interpolation = 1;
    }

    ti = MPI_Wtime() - t0;

    /* fields (norm = volume fraction of the norm phase, vof as set above) and shifts to the database */
    t0 = MPI_Wtime();

    for(i = 0; i < g_norm_count; i++){

        int c = i * g_norm_coarse;

        g_norms[i] = (c < n_int) ? g_vof[c][g_norm_phase] : 0.0;
    }

    err  = rcfd_set_fields(h, g_norms, NULL, NULL);
    err |= rcfd_write_fields(h);
    err |= rcfd_write_shifts(h);

    if(err != 0) of_abort("writing fields/shifts of frame %d failed", i_frame);

    for(p = 0; p < np; p++){

        long long a, b;

        prep_count_written_shifts(p, &a, &b);

        loc[RCFD_OF_MAX_PHASES + p]     = a + b;
        loc[2 * RCFD_OF_MAX_PHASES + p] = b;
    }

    tw = MPI_Wtime() - t0;

    MPI_Reduce(loc, glob, 3 * RCFD_OF_MAX_PHASES, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(dl, dg, np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ml, mg, np, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(of_rank == 0){

        char line[2048];
        int  len = 0;

        len += snprintf(line + len, sizeof(line) - (size_t)len, "[rcfd_of] frame %4d/%d:", i_frame, g_n_frames);

        for(p = 0; p < np && len < (int)sizeof(line) - 200; p++){

            long long prim = glob[p], tot = glob[RCFD_OF_MAX_PHASES + p], inc = glob[2 * RCFD_OF_MAX_PHASES + p];

            len += snprintf(line + len, sizeof(line) - (size_t)len,
                            " phase %d: %lld primary, %lld total shifts (%lld from other ranks), displacement rms/max %.2f/%.2f mm |",
                            p, prim, tot, inc, 1.0e3 * sqrt(dg[p] / (prim > 0 ? (double)prim : 1.0)), 1.0e3 * sqrt(mg[p]));

            g_total_shifts[p] += tot; g_total_incoming[p] += inc;
        }

        printf("%s primary %.3f s, interpolation %.3f s, write %.3f s\n", line, ts, ti, tw);
        fflush(stdout);
    }

    g_t_shift += ts; g_t_interp += ti; g_t_write += tw;

    g_frames_written++;

    return g_frames_written;
}

static int bridge_close(void)
{
    char    name[256];
    int     p, missing = 0, g_missing = 0, complete, vof_files = 1;
    double  t0;

    if(!g_open) return -1;

    complete = (g_frames_written == g_n_frames);

    if(!complete && g_frames_written < 1)
        of_msg("WARNING: no frame was written. There is no database; rerun with the recording window inside the simulated time.");

    /* recurrence path */
    t0 = MPI_Wtime();

    if(complete && rcfd_write_recurrence(g_h) != 0) of_abort("rcfd_write_recurrence failed");

    /* The run ended before all frames were recorded: the frames written are valid. They get a
     * replay path (frame i is followed by i + 1, the last one by the first) so that they can
     * be used; the recurrence format of the header needs all frames and is not applied. */
    if(!complete && g_frames_written >= 1){

        if(of_rank == 0){

            FILE *f;
            int  i;

            bridge_mkdir("data/rec");

            f = fopen("data/rec/jump_0_0_0", "w");

            if(f == NULL) of_abort("cannot write 'data/rec/jump_0_0_0': %s", strerror(errno));

            for(i = 0; i < g_frames_written; i++) fprintf(f, "%d \n", (i + 1) % g_frames_written);

            fclose(f);
        }

        of_msg("WARNING: only %d of %d frames were written (the run ended before the recording was complete). "
               "The database holds these %d frames with a replay path: set number_of_frames to %d in rCFD/rCFD_user_run.h to use it.",
               g_frames_written, g_n_frames, g_frames_written, g_frames_written);

        g_n_frames = g_frames_written;

        bridge_write_prep_info();
    }

    g_t_write += MPI_Wtime() - t0;

    free(g_displ); free(g_norms); free(g_vof); free(g_vof_mem);

    g_displ = NULL; g_norms = NULL; g_vof = NULL; g_vof_mem = NULL;

    vof_files = rcfd_is_vof_on(g_h);        /* a single phase without vof has no vof files */

    /* writes data/rcfd_db.meta.toml */
    rcfd_close(g_h);

    g_h = NULL;

    if(complete){

        for(p = 0; p < g_n_phases; p++){

            snprintf(name, sizeof(name), "data/c2c/c2c_0_%d_%d", p, of_rank);
            if(!of_file_exists(name)) missing++;

            snprintf(name, sizeof(name), "data/vof/vof_0_%d_%d", p, of_rank);
            if(vof_files && !of_file_exists(name)) missing++;
        }

        if(of_rank == 0){

            if(!of_file_exists("data/rec/jump_0_0_0"))     missing++;
            if(!of_file_exists("data/rcfd_db.meta.toml"))  missing++;
        }

        MPI_Allreduce(&missing, &g_missing, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if(g_missing > 0) of_abort("%d database files are missing after the preparation (all ranks)", g_missing);

        if(of_rank == 0){

            printf("[rcfd_of] database written: %d frames;", g_n_frames);

            for(p = 0; p < g_n_phases; p++) printf(" phase %d: %lld shifts (%lld from other ranks);", p, g_total_shifts[p], g_total_incoming[p]);

            printf(" primary shifts %.2f s, interpolation %.2f s, write %.2f s (rank 0)\n", g_t_shift, g_t_interp, g_t_write);
            fflush(stdout);
        }

        of_msg("peak resident memory per rank: max %.0f MB, sum %.0f MB", of_peak_memory_mb(0), of_peak_memory_mb(1));
    }

    of_free_topo(&of_topo);

    /* the user library stays loaded: the library may still hold pointers into it */

    g_open = 0;

    if(g_mpi_started_here){ MPI_Finalize(); g_mpi_started_here = 0; }

    return complete ? 0 : 1;
}

/* The calling program may trap floating point exceptions (OpenFOAM does with
 * FOAM_SIGFPE). The library is not written for that, so the traps are off while it works. */
int rcfd_of_prep_open(const rcfd_of_mesh_t *mesh, const rcfd_of_prep_settings_t *s, rcfd_of_prep_info_t *info)
{
    fenv_t env;
    int    rc;

    feholdexcept(&env);
    rc = bridge_open(mesh, s, info);
    fesetenv(&env);

    return rc;
}

int rcfd_of_prep_frame(const double *vof, const double *origin)
{
    fenv_t env;
    int    rc;

    feholdexcept(&env);
    rc = bridge_frame(vof, origin);
    fesetenv(&env);

    return rc;
}

int rcfd_of_prep_close(void)
{
    fenv_t env;
    int    rc;

    feholdexcept(&env);
    rc = bridge_close();
    fesetenv(&env);

    return rc;
}
