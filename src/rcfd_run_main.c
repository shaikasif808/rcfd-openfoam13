/* rcfd_run - rCFD run mode on a database written during an OpenFOAM run (TUHH 11/2026)
 *
 *   mpirun -np N rCFD/rcfd_run [-steps 800] [-snap 10] [-case DIR]
 *
 * Counterpart of the Fluent run UDF. This program holds no model: phases, number of data
 * (species), drift, diffusion, balances, initial state, sources and monitors come from the
 * user header of the case, rCFD/rCFD_user_run.h, which rcfd_compile builds into this
 * program. Its user_register_functions(handle) is called before the library initialises.
 *
 * Reads   data/topo/{meta.txt,topo_<rank>,gid_<rank>.bin}, data/{c2c,vof,rec}, data/rcfd_db.meta.toml
 * Writes  out/snap_%06d.bin  volume fraction of every phase, every data of every phase and the
 *                            arrays made with rcfd_user_field(), whole mesh, OpenFOAM cell order
 *         out/fields.txt     names of the fields in a snapshot
 *         whatever the user header and the library write (monitors, rCFD_balance.out, ...)
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "rcfd_of_common.h"

#define RUN_MAX_PHASES       8
#define RUN_MAX_USER_FIELDS  64

/* rcfd_user_wrap.c with the user header of the case */
void rcfd_user_register_functions(void *handle);
unsigned long rcfd_user_header_crc(void);
long          rcfd_user_header_bytes(void);

/* ------------------------------------------------------------------------- */
/* boundary patches (data/topo/patches_<rank>, written by rCFDPrep)            */
/* ------------------------------------------------------------------------- */
#define RUN_MAX_PATCHES 256

typedef struct {
    char    name[128], type[64];
    int     n_faces;
    int     *cell;          /* [n_faces] cell next to the face (index into the cell arrays) */
    double  *area;          /* [n_faces] face area [m2] */
    double  *centre;        /* [3 * n_faces] face centre [m] */
} run_patch_t;

static run_patch_t  g_patch[RUN_MAX_PATCHES];
static int          g_n_patches = -1;       /* -1: the database has no patch file */

static void run_read_patches(int n_int)
{
    char    name[256];
    FILE    *f;
    int     i, k, n, ok, g_ok;

    snprintf(name, sizeof(name), "data/topo/patches_%d", of_rank);

    f = fopen(name, "r");

    ok = (f != NULL);

    MPI_Allreduce(&ok, &g_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if(!g_ok){ if(f) fclose(f); g_n_patches = -1; return; }

    if(fscanf(f, "%d", &n) != 1 || n < 0 || n > RUN_MAX_PATCHES) of_abort("'%s': bad number of patches", name);

    for(i = 0; i < n; i++){

        run_patch_t *P = &g_patch[i];

        if(fscanf(f, "%127s %63s %d", P->name, P->type, &P->n_faces) != 3 || P->n_faces < 0)
            of_abort("'%s': bad header of patch %d", name, i);

        P->cell   = (int*)   of_malloc((size_t)P->n_faces * sizeof(int), "patch cells");
        P->area   = (double*)of_malloc((size_t)P->n_faces * sizeof(double), "patch areas");
        P->centre = (double*)of_malloc((size_t)P->n_faces * 3 * sizeof(double), "patch centres");

        for(k = 0; k < P->n_faces; k++){

            if(fscanf(f, "%d %lf %lf %lf %lf", &P->cell[k], &P->area[k], &P->centre[3 * k], &P->centre[3 * k + 1], &P->centre[3 * k + 2]) != 5)
                of_abort("'%s': patch %s, face %d of %d cannot be read", name, P->name, k, P->n_faces);

            if(P->cell[k] < 0 || P->cell[k] >= n_int)
                of_abort("'%s': patch %s, face %d: cell %d is not a cell of this rank (%d cells)", name, P->name, k, P->cell[k], n_int);
        }
    }

    fclose(f);

    g_n_patches = n;
}

/* udf.h: faces of a boundary patch on this rank. Returns the number of faces (0 if this rank
 * has none), -1 if the database has no patch of that name. cells, areas, centres may be NULL. */
int rcfd_run_patch_faces(void *handle, const char *patch, const int **cells, const double **areas, const double **centres)
{
    static int told = 0;
    int i;

    (void)handle;

    if(cells) *cells = NULL;
    if(areas) *areas = NULL;
    if(centres) *centres = NULL;

    if(g_n_patches < 0){

        if(!told) of_msg("rcfd_patch_faces(\"%s\"): the database has no patch information (data/topo/patches_<rank>); it was prepared with an older rcfd_foam, run the preparation again", patch);

        told = 1;

        return -1;
    }

    for(i = 0; i < g_n_patches; i++){

        if(strcmp(g_patch[i].name, patch) == 0){

            if(cells) *cells = g_patch[i].cell;
            if(areas) *areas = g_patch[i].area;
            if(centres) *centres = g_patch[i].centre;

            return g_patch[i].n_faces;
        }
    }

    if(!told){

        char   list[2048] = "";
        size_t used = 0;

        for(i = 0; i < g_n_patches && used + 2 < sizeof(list); i++)
            used += (size_t)snprintf(list + used, sizeof(list) - used, " %.100s", g_patch[i].name);

        of_msg("rcfd_patch_faces(\"%s\"): no such patch. Patches of the database:%s", patch, list);

        told = 1;
    }

    return -1;
}

/* rCFD_user_run.h next to this executable: stop if it was changed after the build */
static void run_check_user_header(const char *argv0)
{
    char    exe[4096], header[4200], *slash;
    ssize_t n = readlink("/proc/self/exe", exe, sizeof(exe) - 1);

    if(n > 0) exe[n] = 0;
    else      snprintf(exe, sizeof(exe), "%s", argv0 ? argv0 : "");

    slash = strrchr(exe, '/');

    if(slash == NULL) return;

    snprintf(header, sizeof(header), "%.*s/rCFD_user_run.h", (int)(slash - exe), exe);

    of_check_user_header(header, exe, rcfd_user_header_crc(), rcfd_user_header_bytes(), "rcfd_compile run");
}

/* ------------------------------------------------------------------------- */
/* user functions                                                             */
/* ------------------------------------------------------------------------- */

static void (*g_user_solver_dict)(void *) = NULL;

static void run_set_solver_dict(void *h)
{
    if(g_user_solver_dict) g_user_solver_dict(h);

    /* No effect on the run itself: with 0 the library allocates the tracer random-walk
       switches of the Prep_Dict with sizeof(short) (rCFD_init.h:559) and clears them with
       sizeof(int) (rCFD_defaults.h:230). With 1 these arrays are not allocated. */
    rcfd_set_dict_property_i(h, "Solver_Dict", -1, -1, -1, "shifts_from_displacement_fields_on", 1);
}

int rcfd_run_register_user_function(void *handle, const char *name, void (*f)(void))
{
    if(strcmp(name, "user_set_solver_dict") == 0){

        g_user_solver_dict = (void(*)(void*))f;

        return 0;       /* run_set_solver_dict is registered and calls it */
    }

    /* run mode: the library reads cells, faces and hosts from data/topo itself */
    if(strcmp(name, "user_set_cell_dict") == 0 || strcmp(name, "user_set_face_dict") == 0 ||
       strcmp(name, "user_set_node_dict") == 0 || strcmp(name, "user_init_cells") == 0 ||
       strcmp(name, "user_init_faces") == 0 || strcmp(name, "user_init_cell_hosts") == 0 ||
       strcmp(name, "user_init_nodes") == 0){

        of_msg("note: %s of the user header is not used in run mode; the topology is read from data/topo", name);

        return 0;
    }

    if(rcfd_register_user_function(handle, name, f) != 0){

        of_msg("WARNING: the library does not know the user function '%s' (ignored)", name);

        return -1;
    }

    return 0;
}

/* named per-cell arrays of the user header (counterpart of C_UDMI) */
static int     g_n_user_fields = 0;
static char    g_user_field_name[RUN_MAX_USER_FIELDS][64];
static double  *g_user_field[RUN_MAX_USER_FIELDS];

double *rcfd_run_user_field(void *handle, const char *name)
{
    int i;

    (void)handle;

    if(name == NULL || name[0] == 0) return NULL;

    for(i = 0; i < g_n_user_fields; i++){ if(strcmp(g_user_field_name[i], name) == 0) return g_user_field[i]; }

    if(g_n_user_fields >= RUN_MAX_USER_FIELDS) of_abort("rcfd_user_field: more than %d fields", RUN_MAX_USER_FIELDS);

    if(strlen(name) > 63) of_abort("rcfd_user_field: name '%s' is longer than 63 characters", name);

    for(i = 0; name[i]; i++){

        if(!(isalnum((unsigned char)name[i]) || name[i] == '_' || name[i] == '.'))
            of_abort("rcfd_user_field: name '%s' may hold letters, digits, '_' and '.' only", name);
    }

    strcpy(g_user_field_name[g_n_user_fields], name);

    g_user_field[g_n_user_fields] = (double*)calloc((size_t)(of_topo.n_int > 0 ? of_topo.n_int : 1), sizeof(double));

    if(g_user_field[g_n_user_fields] == NULL) of_abort("out of memory (user field '%s')", name);

    return g_user_field[g_n_user_fields++];
}

/* ------------------------------------------------------------------------- */
/* database checks before the library reads it                                */
/* ------------------------------------------------------------------------- */

static int run_key_int(const char *text, const char *key)
{
    const char *p = text;
    size_t      len = strlen(key);

    while(p && *p){

        const char *q = p;

        while(*q == ' ' || *q == '\t') q++;

        if(strncmp(q, key, len) == 0){

            q += len;

            while(*q == ' ' || *q == '\t') q++;

            if(*q == '=') return atoi(q + 1);
        }

        p = strchr(p, '\n');
        if(p) p++;
    }

    return -1;
}

static int run_key_text(const char *text, const char *key, char *value, size_t size)
{
    const char *p = text;
    size_t      len = strlen(key);

    value[0] = 0;

    while(p && *p){

        const char *q = p;

        while(*q == ' ' || *q == '\t') q++;

        if(strncmp(q, key, len) == 0){

            q += len;

            while(*q == ' ' || *q == '\t') q++;

            if(*q == '='){

                size_t n = 0;

                q++;

                while(*q == ' ' || *q == '\t') q++;

                while(*q && *q != '\n' && *q != '\r' && n + 1 < size) value[n++] = *q++;

                while(n > 0 && (value[n - 1] == ' ' || value[n - 1] == '\t')) n--;

                value[n] = 0;

                return 1;
            }
        }

        p = strchr(p, '\n');
        if(p) p++;
    }

    return 0;
}

/* ------------------------------------------------------------------------- */
/* snapshots                                                                  */
/* ------------------------------------------------------------------------- */

static int   *g_counts = NULL, *g_displs = NULL;    /* rank 0: internal cells per rank */
static int   *g_gid_all = NULL;                     /* rank 0: global id of every gathered cell */
static float *g_snap_local = NULL;                  /* [n_int] */
static float *g_snap_all = NULL;                    /* rank 0: [n_global] gathered, rank order */
static float *g_snap_field = NULL;                  /* rank 0: [n_global] global order */

static int   g_n_phases = 0, g_n_data[RUN_MAX_PHASES];
static char  g_phase_name[RUN_MAX_PHASES][64];

static void run_output_init(void)
{
    int n_int = of_topo.n_int, r;

    if(of_rank == 0){

        g_counts = (int*)of_malloc((size_t)of_size * sizeof(int), "counts");
        g_displs = (int*)of_malloc((size_t)of_size * sizeof(int), "displs");
    }

    MPI_Gather(&n_int, 1, MPI_INT, g_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(of_rank == 0){

        g_displs[0] = 0;

        for(r = 1; r < of_size; r++) g_displs[r] = g_displs[r - 1] + g_counts[r - 1];

        g_gid_all    = (int*)  of_malloc((size_t)of_n_cells_global * sizeof(int), "gid map");
        g_snap_all   = (float*)of_malloc((size_t)of_n_cells_global * sizeof(float), "snapshot gather");
        g_snap_field = (float*)of_malloc((size_t)of_n_cells_global * sizeof(float), "snapshot field");
    }

    MPI_Gatherv(of_topo.gid, n_int, MPI_INT, g_gid_all, g_counts, g_displs, MPI_INT, 0, MPI_COMM_WORLD);

    g_snap_local = (float*)of_malloc((size_t)n_int * sizeof(float), "snapshot buffer");

    if(of_rank == 0){

        if(mkdir("out", 0777) != 0 && errno != EEXIST) of_abort("cannot create directory 'out': %s", strerror(errno));
    }
}

static int run_field_count(void)
{
    int p, n = 0;

    for(p = 0; p < g_n_phases; p++) n += 1 + g_n_data[p];

    return n + g_n_user_fields;
}

/* out/fields.txt: one name per line, in the order of the snapshot */
static void run_write_field_names(void)
{
    FILE *f;
    int  p, d, i;

    if(of_rank != 0) return;

    f = fopen("out/fields.txt", "w");

    if(f == NULL) of_abort("cannot write 'out/fields.txt': %s", strerror(errno));

    for(p = 0; p < g_n_phases; p++){

        fprintf(f, "vof.%s\n", g_phase_name[p]);

        for(d = 0; d < g_n_data[p]; d++) fprintf(f, "c%d.%s\n", d, g_phase_name[p]);
    }

    for(i = 0; i < g_n_user_fields; i++) fprintf(f, "%s\n", g_user_field_name[i]);

    fclose(f);
}

static void run_gather_and_write(FILE *f, const char *name)
{
    size_t n = (size_t)of_n_cells_global;
    int    i;

    MPI_Gatherv(g_snap_local, of_topo.n_int, MPI_FLOAT, g_snap_all, g_counts, g_displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if(of_rank == 0){

        for(i = 0; i < of_n_cells_global; i++) g_snap_field[g_gid_all[i]] = g_snap_all[i];

        if(of_fwrite_le(g_snap_field, sizeof(float), n, f) != n) of_abort("write error on '%s'", name);
    }
}

/* out/snap_%06d.bin: int64 n_cells, int32 n_fields, float64 time, then n_fields float32 arrays */
static void run_write_snapshot(void *h, int step, double time)
{
    int     i, p, d, k, n_int = of_topo.n_int, i_frame = rcfd_get_current_frame(h);
    double  **vof = rcfd_get_cell_vof(h, 0, i_frame);
    char    name[256], tmp[256];
    FILE    *f = NULL;

    /* a single phase without vof: the library holds no volume fractions, the phase fills the cells */
    if(vof == NULL && rcfd_is_vof_on(h)) of_abort("snapshot: the library returned no volume fractions (frame %d)", i_frame);

    snprintf(name, sizeof(name), "out/snap_%06d.bin", step);
    snprintf(tmp,  sizeof(tmp),  "out/snap_%06d.bin.tmp", step);

    if(of_rank == 0){

        int64_t n64 = of_n_cells_global;
        int32_t nf  = run_field_count();
        double  t   = time;

        f = fopen(tmp, "wb");

        if(f == NULL) of_abort("cannot write '%s': %s", tmp, strerror(errno));

        if(of_fwrite_le(&n64, sizeof(int64_t), 1, f) != 1 || of_fwrite_le(&nf, sizeof(int32_t), 1, f) != 1 || of_fwrite_le(&t, sizeof(double), 1, f) != 1)
            of_abort("write error on '%s'", tmp);
    }

    for(p = 0; p < g_n_phases; p++){

        double **data = rcfd_get_cell_data(h, 0, p);

        for(i = 0; i < n_int; i++) g_snap_local[i] = vof ? (float)vof[i][p] : 1.0f;

        run_gather_and_write(f, tmp);

        for(d = 0; d < g_n_data[p]; d++){

            if(data == NULL) of_abort("snapshot: the library returned no data of phase %d", p);

            for(i = 0; i < n_int; i++) g_snap_local[i] = (float)data[i][d];

            run_gather_and_write(f, tmp);
        }
    }

    for(k = 0; k < g_n_user_fields; k++){

        for(i = 0; i < n_int; i++) g_snap_local[i] = (float)g_user_field[k][i];

        run_gather_and_write(f, tmp);
    }

    if(of_rank == 0){

        if(fclose(f) != 0) of_abort("write error on '%s'", tmp);

        if(rename(tmp, name) != 0) of_abort("cannot rename '%s' to '%s': %s", tmp, name, strerror(errno));
    }

    run_write_field_names();
}

/* ------------------------------------------------------------------------- */

static void run_usage(void)
{
    if(of_rank == 0)
        fprintf(stderr, "usage: mpirun -np N rCFD/rcfd_run [-steps N] [-snap K] [-case DIR]\n"
                        "       -steps N   number of rCFD time steps [100]\n"
                        "       -snap K    snapshot every K steps [10]; 0: first and last only\n"
                        "       -case DIR  case directory [current directory]\n");
}

int main(int argc, char **argv)
{
    void    *h = NULL;
    int     i, p, step, n_steps = 100, snap_every = 10, last_snap = -1;
    int     db_frames = 0, db_partitions = 0, n_frames, prep_phases = -1, prep_frames = -1;
    int     missing = 0, g_missing = 0;
    double  dt, prep_dt = -1.0, t_start, t_run;
    char    *text, phases[512] = "", name[256];
    const char *case_dir = NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &of_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &of_size);

    (void)of_fread_le; (void)of_read_frames_meta; (void)of_read_cfg; (void)of_startup; (void)of_register_topology;
    (void)of_user_set_node_dict; (void)of_user_init_cells; (void)of_user_init_faces; (void)of_user_init_cell_hosts; (void)of_user_init_nodes;
    (void)of_user_set_cell_dict; (void)of_user_set_face_dict;

    for(i = 1; i < argc; i++){

        if(strcmp(argv[i], "-steps") == 0 && i + 1 < argc)      n_steps = atoi(argv[++i]);
        else if(strcmp(argv[i], "-snap") == 0 && i + 1 < argc)  snap_every = atoi(argv[++i]);
        else if(strcmp(argv[i], "-case") == 0 && i + 1 < argc)  case_dir = argv[++i];
        else{ run_usage(); MPI_Finalize(); return 1; }
    }

    if(n_steps < 0 || snap_every < 0){ run_usage(); MPI_Finalize(); return 1; }

    run_check_user_header(argv[0]);

    if(case_dir && chdir(case_dir) != 0) of_abort("cannot change into the case directory '%s': %s", case_dir, strerror(errno));

    t_start = MPI_Wtime();

    /* topology of this rank (the library reads the same files itself) */
    of_read_topo_meta();
    of_read_topo(&of_topo);
    of_check_topo_global(&of_topo);

    run_read_patches(of_topo.n_int);

    /* database */
    text = of_read_text_bcast("data/rcfd_db.meta.toml");

    if(text == NULL) of_abort("database metadata 'data/rcfd_db.meta.toml' not found: run the preparation (foamRun with rCFDPrep) first");

    db_partitions = run_key_int(text, "partitions");
    db_frames     = run_key_int(text, "frames");

    free(text);

    if(db_partitions != of_size)
        of_abort("the database was prepared on %d ranks, this run uses %d: the shifts are partition specific, use %d ranks", db_partitions, of_size, db_partitions);

    if(db_frames < 1) of_abort("data/rcfd_db.meta.toml: no valid 'frames' entry");

    /* what the preparation was run with (written by rCFDPrep) */
    text = of_read_text_bcast("data/rcfd_foam_prep.txt");

    if(text){

        char val[64];

        prep_phases = run_key_int(text, "n_phases");
        prep_frames = run_key_int(text, "n_frames");

        if(run_key_text(text, "frame_interval", val, sizeof(val))) prep_dt = atof(val);

        run_key_text(text, "phases", phases, sizeof(phases));

        free(text);
    }

    /* library in run mode */
    rcfd_open(0, NULL, MPI_COMM_WORLD, &h);

    rcfd_register_user_function(h, "user_set_solver_dict", (void(*)(void))run_set_solver_dict);

    rcfd_user_register_functions(h);        /* -> user_register_functions(handle) of rCFD_user_run.h */

    rcfd_init(h, RCFD_MODE_RUN);

    if(!rcfd_is_initialized(h)) of_abort("rcfd_init failed");

    of_check_library_topology(h, &of_topo);

    n_frames   = rcfd_get_frame_count(h);
    g_n_phases = rcfd_get_phase_count(h);
    dt         = rcfd_get_global_timestep(h);

    if(g_n_phases < 1 || g_n_phases > RUN_MAX_PHASES) of_abort("rCFD_user_run.h: number_of_phases = %d", g_n_phases);

    if(prep_phases > 0 && g_n_phases != prep_phases)
        of_abort("rCFD_user_run.h sets number_of_phases = %d, the database was prepared with %d (%s)", g_n_phases, prep_phases, phases);

    if(n_frames > db_frames)
        of_abort("rCFD_user_run.h sets number_of_frames = %d, the database holds %d", n_frames, db_frames);

    if(n_frames < db_frames)
        of_msg("note: rCFD_user_run.h uses %d of the %d frames of the database", n_frames, db_frames);

    if(prep_dt > 0.0 && fabs(dt - prep_dt) > 1.0e-9 * prep_dt)
        of_abort("rCFD_user_run.h sets global_time_step = %.12g s, the database was prepared with a frame interval of %.12g s", dt, prep_dt);

    (void)prep_frames;

    for(p = 0; p < g_n_phases; p++){

        g_n_data[p] = rcfd_get_phase_data_count(h, p);

        if(g_n_data[p] < 0) g_n_data[p] = 0;

        snprintf(name, sizeof(name), "data/c2c/c2c_0_%d_%d", p, of_rank);
        if(!of_file_exists(name)) missing++;

        snprintf(name, sizeof(name), "data/vof/vof_0_%d_%d", p, of_rank);
        if(rcfd_is_vof_on(h) && !of_file_exists(name)) missing++;
    }

    if(of_rank == 0 && !of_file_exists("data/rec/jump_0_0_0")) missing++;

    MPI_Allreduce(&missing, &g_missing, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(g_missing > 0) of_abort("%d database files are missing (all ranks)", g_missing);

    /* phase names for the field names: from the preparation, else p0, p1, ... */
    {
        char *tok = strtok(phases, " \t");

        for(p = 0; p < g_n_phases; p++){

            if(tok){ snprintf(g_phase_name[p], sizeof(g_phase_name[p]), "%s", tok); tok = strtok(NULL, " \t"); }
            else   { snprintf(g_phase_name[p], sizeof(g_phase_name[p]), "p%d", p); }
        }
    }

    /* the recurrence path must stay inside the frames used */
    if(of_rank == 0){

        FILE *f = fopen("data/rec/jump_0_0_0", "r");
        int  next = -1;

        if(f == NULL) of_abort("cannot open 'data/rec/jump_0_0_0'");

        for(i = 0; i < n_frames; i++){

            if(fscanf(f, "%d", &next) != 1) of_abort("data/rec/jump_0_0_0 has only %d entries, %d frames are used", i, n_frames);

            if(next < 0 || next >= n_frames)
                of_abort("data/rec/jump_0_0_0: frame %d continues with frame %d, but rCFD_user_run.h uses %d frames (the database holds %d): set number_of_frames = %d",
                         i, next, n_frames, db_frames, db_frames);
        }

        fclose(f);
    }

    /* the cell centres the library read must be the ones this program read */
    {
        double **x = rcfd_get_cell_center(h, 0);
        int    bad = 0;

        for(i = 0; i < of_topo.n_cells; i++){

            if(x[i][0] != of_topo.x[3 * i] || x[i][1] != of_topo.x[3 * i + 1] || x[i][2] != of_topo.x[3 * i + 2]) bad++;
        }

        if(bad > 0) of_abort("%d cell centres held by the library differ from data/topo/topo_%d", bad, of_rank);
    }

    rcfd_read_db(h);

    if(of_rank == 0){

        printf("[rcfd_of] run (rCFD_user_run.h): %d ranks, %d cells, %d frames of %g s; phases:", of_size, of_n_cells_global, n_frames, dt);

        for(p = 0; p < g_n_phases; p++) printf(" %s (%d data)", g_phase_name[p], g_n_data[p]);

        printf("; %d steps, snapshot every %d\n", n_steps, snap_every);
        fflush(stdout);
    }

    run_output_init();

    run_write_snapshot(h, 0, 0.0);
    last_snap = 0;

    t_run = MPI_Wtime();

    for(step = 1; step <= n_steps; step++){

        rcfd_run(h);

        if(snap_every > 0 && (step % snap_every) == 0){

            run_write_snapshot(h, step, step * dt);
            last_snap = step;

            of_msg("step %d of %d, t = %g s", step, n_steps, step * dt);
        }
    }

    if(last_snap != n_steps) run_write_snapshot(h, n_steps, n_steps * dt);

    t_run = MPI_Wtime() - t_run;

    of_msg("run finished: %d steps, t = %g s, %d fields per snapshot (out/fields.txt); run time %.2f s (%.4f s per step), total %.2f s; peak resident memory per rank: max %.0f MB, sum %.0f MB",
           n_steps, n_steps * dt, run_field_count(), t_run, (n_steps > 0) ? t_run / n_steps : 0.0, MPI_Wtime() - t_start, of_peak_memory_mb(0), of_peak_memory_mb(1));

    rcfd_close(h);

    free(g_counts); free(g_displs); free(g_gid_all); free(g_snap_local); free(g_snap_all); free(g_snap_field);

    for(i = 0; i < g_n_user_fields; i++) free(g_user_field[i]);

    of_free_topo(&of_topo);

    MPI_Finalize();

    return 0;
}
