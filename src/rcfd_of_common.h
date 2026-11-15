/* rcfd_of_common.h - shared code of rcfd_of_prep and rcfd_of_run
 *
 * rCFD on OpenFOAM/CFDEM data without Fluent (see doc/SPEC.md).
 *
 *   - case.cfg reader (key = value, '#' comments, defaults per SPEC)
 *   - reader of RUN/data/topo/topo_<rank> and gid_<rank>.bin with consistency checks
 *   - user functions that hand the topology to the rCFD library (preparation mode)
 *   - small helpers (abort with message, little-endian binary i/o, meta files)
 *
 * Plain C99 + MPI. Include it in exactly one translation unit per program
 * (everything is static).
 */
#ifndef RCFD_OF_COMMON_H
#define RCFD_OF_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <mpi.h>

#include "rCFD_library.h"

#define OF_PHASE_GAS    0
#define OF_PHASE_SOLID  1
#define OF_DATA_FINE    0
#define OF_DATA_COARSE  1

#define OF_INIT_UNIFORM        0
#define OF_INIT_TRACER_BOTTOM  1

static int of_rank = 0, of_size = 1;

/* ------------------------------------------------------------------------- */
/* messages                                                                   */
/* ------------------------------------------------------------------------- */

static void of_abort(const char *fmt, ...)
{
    va_list ap;
    char    text[1024];

    /* one write per message, so that the messages of several ranks do not interleave */
    va_start(ap, fmt);
    vsnprintf(text, sizeof(text), fmt, ap);
    va_end(ap);

    fflush(stdout);
    fprintf(stderr, "[rcfd_of] ERROR (rank %d): %s\n", of_rank, text);
    fflush(stderr);

    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
}

/* printed by rank 0 only */
static void of_msg(const char *fmt, ...)
{
    va_list ap;

    if(of_rank != 0) return;

    printf("[rcfd_of] ");
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    printf("\n");
    fflush(stdout);
}

static void *of_malloc(size_t n, const char *what)
{
    void *p = malloc(n > 0 ? n : 1);

    if(p == NULL) of_abort("out of memory (%s, %.1f MB)", what, (double)n / 1048576.0);

    return p;
}

static int of_file_exists(const char *name)
{
    struct stat st;

    return (stat(name, &st) == 0);
}

/* POSIX cksum of a file (the number the shell command "cksum" prints); 0: read */
static int of_posix_cksum(const char *name, unsigned long *crc_out, long *bytes_out)
{
    uint32_t        tab[256], crc = 0, c;
    unsigned char   buf[65536];
    size_t          n, k;
    long            bytes = 0, len;
    int             i, j;
    FILE            *f = fopen(name, "rb");

    if(f == NULL) return -1;

    for(i = 0; i < 256; i++){

        c = (uint32_t)i << 24;

        for(j = 0; j < 8; j++) c = (c & 0x80000000u) ? ((c << 1) ^ 0x04C11DB7u) : (c << 1);

        tab[i] = c;
    }

    while((n = fread(buf, 1, sizeof(buf), f)) > 0){

        for(k = 0; k < n; k++) crc = (crc << 8) ^ tab[((crc >> 24) ^ buf[k]) & 0xFFu];

        bytes += (long)n;
    }

    fclose(f);

    for(len = bytes; len > 0; len >>= 8) crc = (crc << 8) ^ tab[((crc >> 24) ^ (uint32_t)(len & 0xFF)) & 0xFFu];

    *crc_out = (unsigned long)(~crc & 0xFFFFFFFFu);
    *bytes_out = bytes;

    return 0;
}

/* Stops if the user header was changed after it was compiled. "built" is what was built from
 * it, crc/bytes are the values rcfd_compile stored in that build (bytes < 0: none stored).
 * A header that is not there is not an error: a build can be used without its source. */
static void of_check_user_header(const char *header, const char *built, unsigned long crc, long bytes, const char *command)
{
    unsigned long   crc_now = 0;
    long            bytes_now = 0;

    if(of_rank != 0 || bytes < 0) return;

    if(of_posix_cksum(header, &crc_now, &bytes_now) != 0){

        of_msg("note: '%s' not found, '%s' is used as it was built", header, built);
        return;
    }

    if(crc_now != crc || bytes_now != bytes)
        of_abort("'%s' was changed after '%s' was built from it, so the change would have no effect. Run '%s' in the case directory, then start again", header, built, command);
}

static long long of_file_size(const char *name)
{
    struct stat st;

    if(stat(name, &st) != 0) return -1;

    return (long long)st.st_size;
}

/* peak resident memory of the ranks in MB (getrusage; kB on Linux): sum != 0: sum over the ranks, else the largest */
static double of_peak_memory_mb(int sum)
{
    struct rusage ru;
    double loc = 0.0, glob = 0.0;

    if(getrusage(RUSAGE_SELF, &ru) == 0) loc = (double)ru.ru_maxrss / 1024.0;

    MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE, sum ? MPI_SUM : MPI_MAX, MPI_COMM_WORLD);

    return glob;
}

/* ------------------------------------------------------------------------- */
/* little-endian binary i/o                                                   */
/* ------------------------------------------------------------------------- */

static int of_host_is_big_endian(void)
{
    const uint16_t one = 1;

    return (*(const unsigned char*)&one == 0);
}

static void of_swap_bytes(void *data, size_t size, size_t count)
{
    unsigned char *p = (unsigned char*)data, t;
    size_t i, j;

    for(i = 0; i < count; i++, p += size){

        for(j = 0; j < size / 2; j++){ t = p[j]; p[j] = p[size - 1 - j]; p[size - 1 - j] = t; }
    }
}

/* read count items of the given size stored little-endian; returns the number of items read */
static size_t of_fread_le(void *data, size_t size, size_t count, FILE *f)
{
    size_t n = fread(data, size, count, f);

    if(of_host_is_big_endian() && size > 1) of_swap_bytes(data, size, n);

    return n;
}

/* write count items little-endian (data is restored on a big-endian host) */
static size_t of_fwrite_le(void *data, size_t size, size_t count, FILE *f)
{
    size_t n;

    if(of_host_is_big_endian() && size > 1) of_swap_bytes(data, size, count);
    n = fwrite(data, size, count, f);
    if(of_host_is_big_endian() && size > 1) of_swap_bytes(data, size, count);

    return n;
}

/* ------------------------------------------------------------------------- */
/* case.cfg                                                                   */
/* ------------------------------------------------------------------------- */

typedef struct {

    /* keys of doc/SPEC.md */
    double  frame_interval;             /* [0.025] s; time step of preparation and run        */
    int     n_frames;                   /* [from frames/meta.txt]                             */
    int     vertical_axis;              /* [1]                                                */
    double  density_solid;              /* [1040]                                             */
    double  d_fine, d_coarse;           /* [0.0018] [0.0033]                                  */
    double  mf_fine, mf_coarse;         /* [0.6] [0.4]                                        */
    double  d32;                        /* [0.0022]                                           */
    double  k_d;                        /* [2.0e5]  u_i = -k_d (d_i - d32)|d_i - d32|         */
    int     n_steps;                    /* [800]                                              */
    int     snap_every;                 /* [10]                                               */
    int     balance_correction_loops;   /* [10]                                               */
    double  face_swap_diffusion;        /* [0.1] solids data (ref/rCFD_user_run_simple.h)     */
    int     interpolation_on;           /* [1]                                                */
    double  vof_cutoff_limit;           /* [1.0e-6]                                           */
    int     init_mode;                  /* [uniform] | tracer_bottom <height>                 */
    double  init_height;
    int     report_priority;            /* [3] = RCFD_REPORT_WARNING                          */

    /* additional optional keys (not in the SPEC; defaults reproduce SPEC / reference) */
    double  density_gas;                /* [1.225]                                            */
    double  face_swap_diffusion_gas;    /* [1.5]  gas data (ref/rCFD_user_run_simple.h)       */
    int     face_swap_diffusion_on;     /* [1]    library default; 0 switches the face swaps off */
    int     data_drifting_on;           /* [-1]   -1: on if k_d != 0, off if k_d == 0         */
    double  vof_cutoff_limit_run;       /* [= vof_cutoff_limit]                               */
    double  next_neigh_cell_skin_factor;/* [1.5]  preparation, as in ref/prep_mpi_test.c      */
    int     prep_bbox_nodes;            /* [1]    preparation: padded bounding box nodes      */
    int     prep_report_priority;       /* [= report_priority]                                */

    int     has_frame_interval, has_n_frames;

} of_cfg_t;

static of_cfg_t of_cfg;

static void of_cfg_defaults(of_cfg_t *c)
{
    c->frame_interval           = 0.025;
    c->n_frames                 = -1;
    c->vertical_axis            = 1;
    c->density_solid            = 1040.0;
    c->d_fine                   = 0.0018;
    c->d_coarse                 = 0.0033;
    c->mf_fine                  = 0.6;
    c->mf_coarse                = 0.4;
    c->d32                      = 0.0022;
    c->k_d                      = 2.0e5;
    c->n_steps                  = 800;
    c->snap_every               = 10;
    c->balance_correction_loops = 10;
    c->face_swap_diffusion      = 0.1;
    c->interpolation_on         = 1;
    c->vof_cutoff_limit         = 1.0e-6;
    c->init_mode                = OF_INIT_UNIFORM;
    c->init_height              = 0.0;
    c->report_priority          = RCFD_REPORT_WARNING;

    c->density_gas                 = 1.225;
    c->face_swap_diffusion_gas     = 1.5;
    c->face_swap_diffusion_on      = 1;
    c->data_drifting_on            = -1;
    c->vof_cutoff_limit_run        = -1.0;
    c->next_neigh_cell_skin_factor = 1.5;
    c->prep_bbox_nodes             = 1;
    c->prep_report_priority        = -1;

    c->has_frame_interval = 0;
    c->has_n_frames       = 0;
}

static char *of_trim(char *s)
{
    char *e;

    while(*s && isspace((unsigned char)*s)) s++;

    e = s + strlen(s);

    while(e > s && isspace((unsigned char)e[-1])) e--;

    *e = 0;

    return s;
}

static double of_cfg_double(const char *key, const char *val, int line)
{
    char *end = NULL;
    double v;

    errno = 0;
    v = strtod(val, &end);

    if(end == val || *of_trim(end) != 0 || errno != 0)
        of_abort("case.cfg line %d: value '%s' of key '%s' is not a number", line, val, key);

    return v;
}

static int of_cfg_int(const char *key, const char *val, int line)
{
    char *end = NULL;
    long v;

    errno = 0;
    v = strtol(val, &end, 10);

    if(end == val || *of_trim(end) != 0 || errno != 0)
        of_abort("case.cfg line %d: value '%s' of key '%s' is not an integer", line, val, key);

    return (int)v;
}

/* parse the text of a case.cfg (modified in place) */
static void of_cfg_parse(of_cfg_t *c, char *text)
{
    char *line = text, *next, *hash, *eq, *key, *val;
    int   i_line = 0;

    while(line && *line){

        next = strchr(line, '\n');
        if(next){ *next = 0; next++; }
        i_line++;

        hash = strchr(line, '#');
        if(hash) *hash = 0;

        key = of_trim(line);

        if(*key){

            eq = strchr(key, '=');

            if(eq == NULL) of_abort("case.cfg line %d: expected 'key = value', found '%s'", i_line, key);

            *eq = 0;
            val = of_trim(eq + 1);
            key = of_trim(key);

            if(*val == 0) of_abort("case.cfg line %d: key '%s' without a value", i_line, key);

#define OF_KEY_D(name)  else if(strcmp(key, #name) == 0){ c->name = of_cfg_double(key, val, i_line); }
#define OF_KEY_I(name)  else if(strcmp(key, #name) == 0){ c->name = of_cfg_int(key, val, i_line); }

            if(strcmp(key, "frame_interval") == 0){ c->frame_interval = of_cfg_double(key, val, i_line); c->has_frame_interval = 1; }
            else if(strcmp(key, "n_frames") == 0){ c->n_frames = of_cfg_int(key, val, i_line); c->has_n_frames = 1; }
            else if(strcmp(key, "init") == 0){

                if(strcmp(val, "uniform") == 0){

                    c->init_mode = OF_INIT_UNIFORM;
                }
                else if(strncmp(val, "tracer_bottom", 13) == 0 && (val[13] == 0 || isspace((unsigned char)val[13]))){

                    char *h = of_trim(val + 13);

                    if(*h == 0) of_abort("case.cfg line %d: 'init = tracer_bottom <height>' needs the height", i_line);

                    c->init_mode   = OF_INIT_TRACER_BOTTOM;
                    c->init_height = of_cfg_double("init = tracer_bottom", h, i_line);
                }
                else{

                    of_abort("case.cfg line %d: init must be 'uniform' or 'tracer_bottom <height>', found '%s'", i_line, val);
                }
            }
            OF_KEY_I(vertical_axis)
            OF_KEY_D(density_solid)
            OF_KEY_D(d_fine)
            OF_KEY_D(d_coarse)
            OF_KEY_D(mf_fine)
            OF_KEY_D(mf_coarse)
            OF_KEY_D(d32)
            OF_KEY_D(k_d)
            OF_KEY_I(n_steps)
            OF_KEY_I(snap_every)
            OF_KEY_I(balance_correction_loops)
            OF_KEY_D(face_swap_diffusion)
            OF_KEY_I(interpolation_on)
            OF_KEY_D(vof_cutoff_limit)
            OF_KEY_I(report_priority)
            OF_KEY_D(density_gas)
            OF_KEY_D(face_swap_diffusion_gas)
            OF_KEY_I(face_swap_diffusion_on)
            OF_KEY_I(data_drifting_on)
            OF_KEY_D(vof_cutoff_limit_run)
            OF_KEY_D(next_neigh_cell_skin_factor)
            OF_KEY_I(prep_bbox_nodes)
            OF_KEY_I(prep_report_priority)
            else{

                if(of_rank == 0){

                    printf("[rcfd_of] WARNING: case.cfg line %d: unknown key '%s' ignored\n", i_line, key);
                    fflush(stdout);
                }
            }

#undef OF_KEY_D
#undef OF_KEY_I
        }

        line = next;
    }
}

static void of_cfg_check(of_cfg_t *c)
{
    if(c->vertical_axis < 0 || c->vertical_axis > 2) of_abort("case.cfg: vertical_axis must be 0, 1 or 2 (is %d)", c->vertical_axis);
    if(!(c->frame_interval > 0.0))                    of_abort("case.cfg: frame_interval must be positive (is %g)", c->frame_interval);
    if(c->has_n_frames && c->n_frames < 1)            of_abort("case.cfg: n_frames must be at least 1 (is %d)", c->n_frames);
    if(c->n_steps < 0)                                of_abort("case.cfg: n_steps must not be negative (is %d)", c->n_steps);
    if(c->balance_correction_loops < 0)               of_abort("case.cfg: balance_correction_loops must not be negative");
    if(!(c->density_solid > 0.0) || !(c->density_gas > 0.0)) of_abort("case.cfg: densities must be positive");
    if(!(c->d_fine > 0.0) || !(c->d_coarse > 0.0) || !(c->d32 > 0.0)) of_abort("case.cfg: d_fine, d_coarse and d32 must be positive");
    if(c->mf_fine < 0.0 || c->mf_coarse < 0.0)        of_abort("case.cfg: mf_fine and mf_coarse must not be negative");
    if(c->vof_cutoff_limit < 0.0)                     of_abort("case.cfg: vof_cutoff_limit must not be negative");
    if(c->report_priority < 0 || c->report_priority >= RCFD_REPORT_MAX) of_abort("case.cfg: report_priority out of range (0..%d)", RCFD_REPORT_MAX - 1);

    if(c->vof_cutoff_limit_run < 0.0) c->vof_cutoff_limit_run = c->vof_cutoff_limit;
    if(c->prep_report_priority < 0)   c->prep_report_priority = c->report_priority;
    if(c->data_drifting_on < 0)       c->data_drifting_on     = (c->k_d != 0.0) ? 1 : 0;
    if(c->snap_every < 0)             c->snap_every           = 0;

    if(fabs(c->mf_fine + c->mf_coarse - 1.0) > 1.0e-9 && of_rank == 0){

        printf("[rcfd_of] WARNING: mf_fine + mf_coarse = %.12g, not 1\n", c->mf_fine + c->mf_coarse);
    }
}

/* rank 0 reads a small text file, all ranks get its content; returns NULL if the file does not exist */
static char *of_read_text_bcast(const char *name)
{
    long  n = -1;
    char *text = NULL;

    if(of_rank == 0){

        FILE *f = fopen(name, "rb");

        if(f){

            if(fseek(f, 0, SEEK_END) == 0) n = ftell(f);

            if(n >= 0){

                rewind(f);
                text = (char*)of_malloc((size_t)n + 1, name);

                if(fread(text, 1, (size_t)n, f) != (size_t)n) of_abort("failed to read '%s'", name);

                text[n] = 0;
            }

            fclose(f);
        }
    }

    MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    if(n < 0) return NULL;

    if(of_rank != 0) text = (char*)of_malloc((size_t)n + 1, name);

    MPI_Bcast(text, (int)n + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

    return text;
}

/* read RUN/case.cfg (current directory must be RUN) */
static void of_read_cfg(void)
{
    char *text = of_read_text_bcast("case.cfg");

    of_cfg_defaults(&of_cfg);

    if(text == NULL) of_abort("cannot open 'case.cfg' in the run directory");

    of_cfg_parse(&of_cfg, text);
    of_cfg_check(&of_cfg);

    free(text);
}

/* ------------------------------------------------------------------------- */
/* meta files                                                                 */
/* ------------------------------------------------------------------------- */

static int of_n_cells_global = 0;       /* data/topo/meta.txt */
static int of_meta_ranks     = 0;

/* data/topo/meta.txt: "n_cells_global n_ranks" */
static void of_read_topo_meta(void)
{
    char *text = of_read_text_bcast("data/topo/meta.txt");
    long long n = 0;
    int r = 0;

    if(text == NULL) of_abort("cannot open 'data/topo/meta.txt'");

    if(sscanf(text, "%lld %d", &n, &r) != 2) of_abort("data/topo/meta.txt: expected 'n_cells_global n_ranks'");

    free(text);

    if(n < 1 || n > 2147483647LL) of_abort("data/topo/meta.txt: n_cells_global = %lld out of range", n);

    if(r != of_size)
        of_abort("topology was written for %d ranks (data/topo/meta.txt), but this program runs on %d ranks", r, of_size);

    of_n_cells_global = (int)n;
    of_meta_ranks     = r;
}

/* frames/meta.txt: "n_frames n_cells_global frame_interval_seconds"; returns 0 if the file does not exist */
static int of_read_frames_meta(int *n_frames, double *interval)
{
    char *text = of_read_text_bcast("frames/meta.txt");
    long long n = 0;
    int nf = 0;
    double dt = 0.0;

    if(text == NULL) return 0;

    if(sscanf(text, "%d %lld %lf", &nf, &n, &dt) != 3) of_abort("frames/meta.txt: expected 'n_frames n_cells_global frame_interval_seconds'");

    free(text);

    if(nf < 1) of_abort("frames/meta.txt: n_frames = %d", nf);

    if(n != (long long)of_n_cells_global)
        of_abort("frames/meta.txt has %lld cells, data/topo/meta.txt has %d", n, of_n_cells_global);

    *n_frames = nf;
    *interval = dt;

    return 1;
}

/* ------------------------------------------------------------------------- */
/* topology of this rank: data/topo/topo_<rank> and gid_<rank>.bin            */
/*                                                                            */
/* topo_<rank> is the library's own ASCII format (rCFD_topo.h: topo_read):    */
/*   line 1 : n_cells n_faces n_nodes                                         */
/*            n_cells = internal + corona cells, n_faces = internal +         */
/*            partition faces of this rank                                    */
/*   n_cells lines : x y z volume host   (internal cells first: host == rank) */
/*   n_faces lines : c0 c1 Ax Ay Az      (faces between internal cells first; */
/*                   partition faces: c0 internal, c1 corona)                 */
/*   node section  : only if n_nodes > 0; not read here (the library reads it */
/*                   in run mode for its VTK dump)                            */
/* ------------------------------------------------------------------------- */

typedef struct {

    int     n_cells, n_int, n_ext;
    int     n_faces, n_fint, n_fext;
    int     n_nodes;                /* as given in the header of the file */

    double  *x;                     /* [3 * n_cells] */
    double  *vol;                   /* [n_cells]     */
    int     *host;                  /* [n_cells]     */

    int     *c0, *c1;               /* [n_faces]     */
    double  *area;                  /* [3 * n_faces] */

    int     *gid;                   /* [n_cells] global (OpenFOAM) cell index */

} of_topo_t;

static of_topo_t of_topo;

static void of_free_topo(of_topo_t *T)
{
    free(T->x); free(T->vol); free(T->host); free(T->c0); free(T->c1); free(T->area); free(T->gid);
    memset(T, 0, sizeof(*T));
}

static void of_read_topo(of_topo_t *T)
{
    char    name[256];
    FILE    *f;
    int     i, n_ext_started = 0, n_fext_started = 0;
    long long size;

    memset(T, 0, sizeof(*T));

    /* topo_<rank> */
    snprintf(name, sizeof(name), "data/topo/topo_%d", of_rank);

    f = fopen(name, "r");

    if(f == NULL) of_abort("cannot open topology file '%s'", name);

    if(fscanf(f, "%d %d %d", &T->n_cells, &T->n_faces, &T->n_nodes) != 3)
        of_abort("%s: header line must be 'n_cells n_faces n_nodes'", name);

    if(T->n_cells < 1)  of_abort("%s: n_cells = %d (a rank without cells is not supported by the library)", name, T->n_cells);
    if(T->n_faces < 1)  of_abort("%s: n_faces = %d (a rank without faces is not supported by the library)", name, T->n_faces);
    if(T->n_nodes < 0)  of_abort("%s: n_nodes = %d", name, T->n_nodes);

    T->x    = (double*)of_malloc(3 * (size_t)T->n_cells * sizeof(double), "cell centres");
    T->vol  = (double*)of_malloc((size_t)T->n_cells * sizeof(double), "cell volumes");
    T->host = (int*)   of_malloc((size_t)T->n_cells * sizeof(int), "cell hosts");

    for(i = 0; i < T->n_cells; i++){

        if(fscanf(f, "%lf %lf %lf %lf %d", &T->x[3 * i], &T->x[3 * i + 1], &T->x[3 * i + 2], &T->vol[i], &T->host[i]) != 5)
            of_abort("%s: cannot read cell line %d of %d (expected 'x y z volume host')", name, i + 1, T->n_cells);

        if(T->host[i] < 0 || T->host[i] >= of_size)
            of_abort("%s: cell %d has host %d, but there are %d ranks", name, i, T->host[i], of_size);

        if(!(T->vol[i] > 0.0)) of_abort("%s: cell %d has volume %g", name, i, T->vol[i]);

        if(T->host[i] == of_rank){

            if(n_ext_started) of_abort("%s: internal cell %d after a corona cell (internal cells must come first)", name, i);

            T->n_int++;
        }
        else{

            n_ext_started = 1;
        }
    }

    T->n_ext = T->n_cells - T->n_int;

    if(T->n_int < 1) of_abort("%s: no internal cell (host == %d)", name, of_rank);

    T->c0   = (int*)   of_malloc((size_t)T->n_faces * sizeof(int), "face c0");
    T->c1   = (int*)   of_malloc((size_t)T->n_faces * sizeof(int), "face c1");
    T->area = (double*)of_malloc(3 * (size_t)T->n_faces * sizeof(double), "face areas");

    for(i = 0; i < T->n_faces; i++){

        if(fscanf(f, "%d %d %lf %lf %lf", &T->c0[i], &T->c1[i], &T->area[3 * i], &T->area[3 * i + 1], &T->area[3 * i + 2]) != 5)
            of_abort("%s: cannot read face line %d of %d (expected 'c0 c1 Ax Ay Az')", name, i + 1, T->n_faces);

        if(T->c0[i] < 0 || T->c0[i] >= T->n_int)
            of_abort("%s: face %d: c0 = %d is not an internal cell (0..%d)", name, i, T->c0[i], T->n_int - 1);

        if(T->c1[i] < 0 || T->c1[i] >= T->n_cells || T->c1[i] == T->c0[i])
            of_abort("%s: face %d: c1 = %d invalid (c0 = %d, %d cells)", name, i, T->c1[i], T->c0[i], T->n_cells);

        if(T->c1[i] < T->n_int){

            if(n_fext_started) of_abort("%s: face %d between internal cells after a partition face (internal faces must come first)", name, i);

            T->n_fint++;
        }
        else{

            n_fext_started = 1;
        }
    }

    T->n_fext = T->n_faces - T->n_fint;

    fclose(f);

    /* gid_<rank>.bin */
    snprintf(name, sizeof(name), "data/topo/gid_%d.bin", of_rank);

    size = of_file_size(name);

    if(size < 0) of_abort("cannot open '%s'", name);

    if(size != 4LL * T->n_cells)
        of_abort("%s has %lld bytes = %lld cells, but topo_%d has %d cells (internal %d + corona %d)",
                 name, size, size / 4, of_rank, T->n_cells, T->n_int, T->n_ext);

    T->gid = (int*)of_malloc((size_t)T->n_cells * sizeof(int), "global cell ids");

    f = fopen(name, "rb");

    if(f == NULL) of_abort("cannot open '%s'", name);

    {
        int32_t *buf = (int32_t*)of_malloc((size_t)T->n_cells * sizeof(int32_t), "gid buffer");

        if(of_fread_le(buf, sizeof(int32_t), (size_t)T->n_cells, f) != (size_t)T->n_cells) of_abort("failed to read '%s'", name);

        for(i = 0; i < T->n_cells; i++){

            if(buf[i] < 0 || buf[i] >= of_n_cells_global)
                of_abort("%s: entry %d is %d, outside 0..%d (data/topo/meta.txt)", name, i, (int)buf[i], of_n_cells_global - 1);

            T->gid[i] = (int)buf[i];
        }

        free(buf);
    }

    fclose(f);
}

/* every global cell must be an internal cell of exactly one rank; corona cells must agree with their owner */
/* image cells behind periodic patches: the last of_n_images of the n_int cells of this
 * rank (preparation only, see rcfd_of_images.h); they copy cells of the mesh and are not
 * counted here */
static int of_n_images = 0;

static void of_check_topo_global(const of_topo_t *T)
{
    int     *own_local, *own_global, i, bad = 0, g_bad = 0;
    long long n_int_sum = 0, n_int = T->n_int - of_n_images;

    MPI_Allreduce(&n_int, &n_int_sum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if(n_int_sum != (long long)of_n_cells_global)
        of_abort("the internal cells of all ranks add up to %lld, data/topo/meta.txt says %d cells", n_int_sum, of_n_cells_global);

    /* owner rank + 1 of every global cell (0 = nobody, sum over ranks detects double ownership) */
    own_local  = (int*)calloc((size_t)of_n_cells_global, sizeof(int));
    own_global = (int*)calloc((size_t)of_n_cells_global, sizeof(int));

    if(own_local == NULL || own_global == NULL) of_abort("out of memory (ownership check)");

    for(i = 0; i < T->n_int - of_n_images; i++){

        if(own_local[T->gid[i]] != 0) of_abort("gid_%d.bin: global cell %d appears twice among the internal cells", of_rank, T->gid[i]);

        own_local[T->gid[i]] = of_rank + 1;
    }

    MPI_Allreduce(own_local, own_global, of_n_cells_global, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /* with n_int_sum == n_global, a cell owned twice implies a cell owned by nobody */
    for(i = 0; i < of_n_cells_global; i++){ if(own_global[i] == 0) bad++; }

    for(i = T->n_int; i < T->n_cells; i++){

        if(own_global[T->gid[i]] != T->host[i] + 1) bad++;
    }

    MPI_Allreduce(&bad, &g_bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    free(own_local);
    free(own_global);

    if(g_bad > 0)
        of_abort("topology/gid files are inconsistent: %d global cells without owner or corona cells whose host does not own the cell (all ranks)", g_bad);
}

/* ------------------------------------------------------------------------- */
/* user functions that hand the topology to the library                       */
/*                                                                            */
/* Preparation mode: the library does not read data/topo itself, all six      */
/* functions below are needed.                                                */
/* Run mode: the library reads data/topo/topo_<rank> itself (init_topo_dict,  */
/* init_topo -> topo_read, parallel_read_hosts); the functions are NOT        */
/* registered there, see of_register_topology().                              */
/* ------------------------------------------------------------------------- */

static double of_bbox_pad = 0.0;

static void of_user_set_cell_dict(void *h, const short i_layer)
{
    if(i_layer != 0) return;

    rcfd_set_dict_property_i(h, "Cell_Dict", i_layer, -1, -1, "number_of_int_cells", of_topo.n_int);
    rcfd_set_dict_property_i(h, "Cell_Dict", i_layer, -1, -1, "number_of_ext_cells", of_topo.n_ext);
    rcfd_set_dict_property_r(h, "Cell_Dict", i_layer, -1, -1, "next_neigh_cell_skin_factor", of_cfg.next_neigh_cell_skin_factor);
}

static void of_user_set_face_dict(void *h, const short i_layer)
{
    if(i_layer != 0) return;

    rcfd_set_dict_property_i(h, "Face_Dict", i_layer, -1, -1, "number_of_int_faces", of_topo.n_fint);
    rcfd_set_dict_property_i(h, "Face_Dict", i_layer, -1, -1, "number_of_ext_faces", of_topo.n_fext);
}

/* Node data are not needed to prepare the database. With number_of_nodes == 0 the
 * library takes the bounding box of a rank from the cell CENTRES (topo_init_bounding_box)
 * and tests the origin points against it with strict inequalities when it looks for
 * primary shifts on other ranks (prep_shift_from_displacement_fields). An origin point
 * within half a cell of a wall is then outside the box of the rank that owns the cell it
 * lies in, and the shift falls back to the nearest cell of the destination rank.
 * Eight pseudo nodes, the corners of the centre box padded by the largest cell size of
 * the rank, give the library a box that contains all cells of the rank. */
static void of_user_set_node_dict(void *h, const short i_layer)
{
    if(i_layer != 0) return;

    rcfd_set_dict_property_i(h, "Node_Dict", i_layer, -1, -1, "number_of_nodes", of_cfg.prep_bbox_nodes ? 8 : 0);
}

static void of_user_init_cells(void *h, const short i_layer)
{
    int     i;
    double  **x  = rcfd_get_cell_center(h, i_layer);
    double  *vol = rcfd_get_cell_volume(h, i_layer);

    if(i_layer != 0) return;

    if(x == NULL || vol == NULL) of_abort("user_init_cells: the library has no cell arrays");

    for(i = 0; i < of_topo.n_cells; i++){

        x[i][0] = of_topo.x[3 * i];
        x[i][1] = of_topo.x[3 * i + 1];
        x[i][2] = of_topo.x[3 * i + 2];
        vol[i]  = of_topo.vol[i];
    }
}

static void of_user_init_faces(void *h, const short i_layer)
{
    int     i;
    int     *c0 = rcfd_get_face_c0(h, i_layer), *c1 = rcfd_get_face_c1(h, i_layer);
    double  **A = rcfd_get_face_area(h, i_layer);

    if(i_layer != 0) return;

    if(c0 == NULL || c1 == NULL || A == NULL) of_abort("user_init_faces: the library has no face arrays");

    for(i = 0; i < of_topo.n_faces; i++){

        c0[i]   = of_topo.c0[i];
        c1[i]   = of_topo.c1[i];
        A[i][0] = of_topo.area[3 * i];
        A[i][1] = of_topo.area[3 * i + 1];
        A[i][2] = of_topo.area[3 * i + 2];
    }
}

static void of_user_init_cell_hosts(void *h, const short i_layer)
{
    int i, *host = rcfd_get_cell_host(h, i_layer);

    if(i_layer != 0) return;

    if(host == NULL) of_abort("user_init_cell_hosts: the library has no host array");

    for(i = 0; i < of_topo.n_cells; i++) host[i] = of_topo.host[i];
}

static void of_user_init_nodes(void *h, const short i_layer)
{
    int     i, k, d;
    double  lo[3] = {1.0e300, 1.0e300, 1.0e300}, hi[3] = {-1.0e300, -1.0e300, -1.0e300}, pad = 0.0, s;
    double  **xn;

    if(i_layer != 0 || !of_cfg.prep_bbox_nodes) return;

    xn = rcfd_get_node(h, i_layer);

    if(xn == NULL) of_abort("user_init_nodes: the library has no node array");

    for(i = 0; i < of_topo.n_cells; i++){

        for(d = 0; d < 3; d++){

            if(of_topo.x[3 * i + d] < lo[d]) lo[d] = of_topo.x[3 * i + d];
            if(of_topo.x[3 * i + d] > hi[d]) hi[d] = of_topo.x[3 * i + d];
        }

        s = cbrt(of_topo.vol[i]);

        if(s > pad) pad = s;
    }

    of_bbox_pad = pad;

    for(k = 0; k < 8; k++){

        for(d = 0; d < 3; d++) xn[k][d] = ((k >> d) & 1) ? hi[d] + pad : lo[d] - pad;
    }
}

/* prep != 0: preparation mode, the topology comes from the functions above.
 * prep == 0: run mode, nothing is registered: init_topo_dict() takes the counts from the
 *            header of topo_<rank>, init_topo() -> topo_read() reads cells, faces and the
 *            optional node section, parallel_set_host_of_cell() -> parallel_read_hosts()
 *            reads the hosts. Registering the functions as well would only overwrite
 *            these values with the same numbers; the node dict function must not be
 *            registered in run mode at all, because a node count that differs from the
 *            file makes topo_read_nodes() fail. */
static void of_register_topology(void *h, int prep)
{
    if(!prep) return;

    rcfd_register_user_function(h, "user_set_cell_dict",   (void(*)(void))of_user_set_cell_dict);
    rcfd_register_user_function(h, "user_set_face_dict",   (void(*)(void))of_user_set_face_dict);
    rcfd_register_user_function(h, "user_set_node_dict",   (void(*)(void))of_user_set_node_dict);
    rcfd_register_user_function(h, "user_init_cells",      (void(*)(void))of_user_init_cells);
    rcfd_register_user_function(h, "user_init_faces",      (void(*)(void))of_user_init_faces);
    rcfd_register_user_function(h, "user_init_cell_hosts", (void(*)(void))of_user_init_cell_hosts);
    rcfd_register_user_function(h, "user_init_nodes",      (void(*)(void))of_user_init_nodes);
}

/* the library must have ended up with the topology of the files */
static void of_check_library_topology(void *h, const of_topo_t *T)
{
    int nc = rcfd_get_cell_count(h, 0), ni = rcfd_get_int_cell_count(h, 0);
    int nf = rcfd_get_face_count(h, 0), nfi = rcfd_get_int_face_count(h, 0);

    if(nc != T->n_cells || ni != T->n_int || nf != T->n_faces || nfi != T->n_fint)
        of_abort("the library holds %d cells (%d internal), %d faces (%d internal); topo_%d has %d (%d), %d (%d)",
                 nc, ni, nf, nfi, of_rank, T->n_cells, T->n_int, T->n_faces, T->n_fint);
}

/* ------------------------------------------------------------------------- */
/* start-up shared by both programs                                           */
/* ------------------------------------------------------------------------- */

/* MPI is initialised; changes into the run directory, reads case.cfg, the meta files,
 * the topology of this rank and checks it */
static void of_startup(int argc, char **argv, const char *prog)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &of_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &of_size);

    (void)of_fwrite_le; (void)of_fread_le; (void)of_read_frames_meta; (void)of_check_library_topology;   /* not used by every program */

    if(argc < 2){

        if(of_rank == 0) fprintf(stderr, "usage: mpirun -np N %s RUN_DIR\n", prog);

        MPI_Finalize();
        exit(1);
    }

    if(chdir(argv[1]) != 0) of_abort("cannot change into run directory '%s': %s", argv[1], strerror(errno));

    of_read_cfg();
    of_read_topo_meta();
    of_read_topo(&of_topo);
    of_check_topo_global(&of_topo);

    {
        long long loc[4], sum[4];

        loc[0] = of_topo.n_int; loc[1] = of_topo.n_ext; loc[2] = of_topo.n_fint; loc[3] = of_topo.n_fext;

        MPI_Allreduce(loc, sum, 4, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        of_msg("%s: %d ranks, %d cells; all ranks: %lld corona cells, %lld internal faces, %lld partition faces (counted on both sides)",
               prog, of_size, of_n_cells_global, sum[1], sum[2], sum[3]);
    }
}

#endif /* RCFD_OF_COMMON_H */
