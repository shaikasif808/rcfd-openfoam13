/******************************************************************************
 * @internal
 * @file        rCFD_library.c
 * @copyright   (C) 2021
 *              Particulate Flow Modelling,
 *              Johannes Kepler University, Linz, Austria
 * @website     www.particulate-flow.at
 * @authors     Stefan Pirker
 * @authors     Daniel Queteschiner
 *****************************************************************************/

#include "rCFD.h"

#include <stdio.h>
#include <string.h> // memset, strcmp
#include <time.h>   // clock

#include "rCFD_library.h"
#include "rCFD_types.h"
#include "rCFD_memory.h"
#include "rCFD_macros.h"
#include "rCFD_report.h"
#include "rCFD_parallel.h"
#include "rCFD_dict.h"
#include "rCFD_dict_read_toml.h"
#include "rCFD_dict_write_toml.h"
#include "rCFD_dump.h"
#include "rCFD_defaults.h"
#include "rCFD_init.h"
#include "rCFD_layer.h"
#include "rCFD_topo.h"
#include "rCFD_run.h"
#include "rCFD_shift.h"
#include "rCFD_prepare.h"
#include "rCFD_swap.h"
#include "rCFD_balance.h"
#include "rCFD_post.h"
#include "rCFD_free.h"
#include "rCFD_user_defaults.h"


/**
 * @brief get rCFD library version as a string of the form "rCFD major.minor.patch"
 *
 * @param major get major version
 * @param minor get minor version
 * @param patch get patch level
 * @return rCFD library version string
 */
const char* rcfd_version(int *major, int *minor, int *patch)
{
    if(major) *major = RCFD_VERSION_MAJOR;
    if(minor) *minor = RCFD_VERSION_MINOR;
    if(patch) *patch = RCFD_VERSION_PATCH;

    return RCFD_VERSION_STRING;
}

/**
 * @brief parse command line args for configuration filename (-c filename)
 * @param argc  number of command-line arguments
 * @param argv  list of command-line arguments
 * @return configuration filename or NULL if config option is not found
 */
char* get_config_filename_from_args(int argc, char **argv)
{
    int count;

    for(count = 0; count < argc-1; ++count){

        if(strcmp(argv[count], "--") == 0)
            break;

        if(    (strcmp(argv[count], "-c") == 0)
            || (strcmp(argv[count], "/c") == 0)
            || (strcmp(argv[count], "-config") == 0)
            || (strcmp(argv[count], "/config") == 0)
            || (strcmp(argv[count], "--config") == 0)){

            return argv[count+1];
        }
    }

    return NULL;
}

/**
 * @brief create an instance of rCFD and return a pointer to it;
 *        pass in command-line args and MPI communicator to run on
 *
 * @param argc  number of command-line arguments
 * @param argv  list of command-line arguments
 * @param comm  MPI communicator
 * @param ptr   pointer to a void pointer variable which serves as a handle
 */
void rcfd_open(int argc, char **argv, MPI_Comm comm, void **ptr)
{
    rCFD_type *rcfd = init_rcfd(comm, get_config_filename_from_args(argc, argv));

    if(ptr) *ptr = (void*)rcfd;
}

/**
 * @brief create an instance of rCFD and return pointer to it;
 *        caller doesn't know MPI communicator, so use MPI_COMM_WORLD;
 *        intialize MPI if needed
 *
 * @param argc  number of command-line arguments
 * @param argv  list of command-line arguments
 * @param ptr   pointer to a void pointer variable which serves as a handle
 */
void rcfd_open_no_mpi(int argc, char **argv, void **ptr)
{
#if defined(RCFD_LIB_MPI)
    int flag = 0;

    MPI_Initialized(&flag);

    if(!flag){

        int argc = 0;
        char **argv = NULL;

        MPI_Init(&argc,&argv);
    }
#endif

    rcfd_open(argc, argv, MPI_COMM_WORLD, ptr);

#if defined(RCFD_LIB_MPI)
    if(!flag){

        rCFD_type *rcfd = (rCFD_type*)(*ptr);

        rcfd->Parallel->finalize = 1;
    }
#endif
}

/**
 * @brief initialize all data structures for specified mode
 *
 * @param handle    pointer to previously created rCFD instance
 * @param mode      preparation (0) or run (1) mode
 */
void rcfd_init(void *handle, int mode)
{
    rCFD_type *rcfd = (rCFD_type*)handle;

    report_start_timer(rcfd, RCFD_TIMER_GLOBAL);

    switch(mode){

    case RCFD_MODE_PREPARATION:
        init_all(rcfd, mode);
        break;
    case RCFD_MODE_RUN:
        report_start_timer(rcfd, RCFD_TIMER_INIT);
        init_all(rcfd, mode);
        post_init_all(rcfd);
        break;
    }
}

/**
 * @brief get initialization status of rCFD instance
 *
 * @param handle    pointer to previously created rCFD instance
 * @return initialization status of rCFD instance
 */
int rcfd_is_initialized(void *handle)
{
    rCFD_type *rcfd = (rCFD_type*)handle;

    return rcfd->initialized;
}

/**
 * @brief read rCFD database (cell-to-cell shift pattern)
 *
 * @param handle    pointer to previously created rCFD instance
 */
void rcfd_read_db(void *handle)
{
    rCFD_type *rcfd = (rCFD_type*)handle;

    report_start_timer(rcfd, RCFD_TIMER_READ_DB);

    report_start_timer(rcfd, RCFD_TIMER_READ_SHIFT);
    shift_read(rcfd);
    post_read_shift(rcfd);

#if MULTI_LAYER_IMPLEMENTATION_COMPLETE
    shift_create_shifts_for_upper_layers(rcfd);
#endif

    parallel_shift_init(rcfd);

    post_read(rcfd);
}

/* ------------------------------------------------------------------------------
 * TUHH 11/2026: where the cell shifts take their data from (diagnostics, no effect
 * on data).
 *
 * In the first RCFD_C2C_DIAG_STEPS steps, right before the cell shifts of a phase
 * with two or more concentration data, the shift lists of the current frame are
 * walked. Every shift is put into one of three classes by its source cell:
 *   cells of this node (interior shift, source is an internal cell),
 *   corona cells       (interior shift, source index is a corona cell),
 *   other nodes        (incoming shift, source data fetched from the other node).
 * Per node and class: number of shifts, share of the node's solids mass that is
 * filled from the class, concentration of the last concentration data pulled from
 * it, mean vertical offset and largest distance between source and destination
 * cell centre, and the weight share and pulled concentration by distance.
 * After the shift the result is compared with the value expected from the lists.
 * For the steps RCFD_C2C_DIAG_FILE_FIRST .. _LAST the nodes selected by
 * RCFD_C2C_DIAG_FILE_NODE write every shift to rCFD_c2c_step<S>_node<N>.txt.
 * Set RCFD_C2C_DIAG_STEPS to 0 to switch all of it off.
 * ------------------------------------------------------------------------------ */
#ifndef RCFD_C2C_DIAG_STEPS
#define RCFD_C2C_DIAG_STEPS 4
#endif
#ifndef RCFD_C2C_DIAG_FILE_FIRST
#define RCFD_C2C_DIAG_FILE_FIRST 1
#endif
#ifndef RCFD_C2C_DIAG_FILE_LAST
#define RCFD_C2C_DIAG_FILE_LAST 2
#endif
#ifndef RCFD_C2C_DIAG_FILE_NODE
#define RCFD_C2C_DIAG_FILE_NODE(n) (((n) == 1) || ((n) == 3) || ((n) == 5))
#endif
#define RCFD_C2C_DIAG_K     6
#define RCFD_C2C_DIAG_NBIN  5
#define RCFD_C2C_DIAG_NCL   (2 * RCFD_C2C_DIAG_NBIN + 2)
#define RCFD_C2C_DIAG_NSTAT (6 + 3 * RCFD_C2C_DIAG_NCL)
#define RCFD_C2C_DIAG_MAX_NODES 128

static double *rcfd_c2c_pred = NULL;
static int     rcfd_c2c_pred_size = 0, rcfd_c2c_pred_valid = 0, rcfd_c2c_pred_data = 0;

static int rcfd_c2c_diag_active(rCFD_type *rcfd, int i_phase)
{
#if RP_NODE
    if(rcfd->Solver->global_run_counter >= RCFD_C2C_DIAG_STEPS) return 0;
    if(rcfd->Solver->current_layer != 0) return 0;
    if(rcfd->Phase_Dict[i_phase].number_of_concentration_data < 2) return 0;
    if(compute_node_count > RCFD_C2C_DIAG_MAX_NODES) return 0;
    return 1;
#else
    return 0;
#endif
}

static void rcfd_c2c_diag(rCFD_type *rcfd, int i_phase)
{
#if RP_NODE
    static const double bin_edge[RCFD_C2C_DIAG_NBIN - 1] = {0.005, 0.010, 0.020, 0.050};
    static const char  *class_name[3] = {"cells of this node", "corona cells      ", "other nodes       "};

    const int K = RCFD_C2C_DIAG_K;
    int     step, i_layer = 0, i_state, i_frame, i_prev, n_int, n_data, i_dc = -1, i_data, n_conc = 0;
    int     i, k, b, j, i_node, i_dim, c0, c1, wchar, fmt;
    int     n_shift, n_in, n_out, n_send, n_total, n_tr, n_mine = 0, i_global;
    int     *recvcounts = NULL, *displs = NULL, *to_node_count, *out2in;
    int     write_file;
    double  *out_vals = NULL, *g_out = NULL, *g_in = NULL, *in_vals = NULL, *wsum = NULL, *pred = NULL;
    double  st[RCFD_C2C_DIAG_NSTAT], md[3] = {0.0, 0.0, 0.0}, g[RCFD_C2C_DIAG_NSTAT], gmd[3];
    double  w, m, f, d, dx, c_src, vof_prev_src, x0[3], rho;
    char    name[96], line[1024];
    int     len;
    FILE    *fo = NULL;
    C2C_Shift_type *S, *Sin, *Sout;

    rcfd_c2c_pred_valid = 0;

    if(!rcfd_c2c_diag_active(rcfd, i_phase)) return;

    step    = rcfd->Solver->global_run_counter;
    i_state = rcfd->Solver->current_state;
    i_frame = rcfd->Solver->current_frame;
    i_prev  = rcfd->Solver->prev_frame;
    n_int   = _CELL_DICT(i_layer).number_of_int_cells;
    n_data  = rcfd->Phase_Dict[i_phase].number_of_data;
    rho     = rcfd->Phase_Dict[i_phase].density;

    if((i_prev < 0) || (i_prev >= rcfd->Solver_Dict->number_of_frames)) i_prev = i_frame;

    /* the last concentration data of the phase (the coarse fraction in a two-fraction case) */
    for(i_data = 0; i_data < n_data; i_data++){

        if(rcfd->Data_Dict[i_phase][i_data].type == RCFD_DATA_CONCENTRATION){ i_dc = i_data; n_conc++; }
    }
    if(i_dc < 0) return;

    S    = &rcfd->C2C[i_state][i_phase][i_frame][i_layer].shifts;
    Sin  = &rcfd->C2C[i_state][i_phase][i_frame][i_layer].shifts_in;
    Sout = &rcfd->C2C[i_state][i_phase][i_frame][i_layer].shifts_out;
    to_node_count = rcfd->C2C[i_state][i_phase][i_frame][i_layer].shifts_to_node_count;
    out2in        = rcfd->C2C[i_state][i_phase][i_frame][i_layer].out2in;

    fmt   = rcfd->Solver->c2c_shift_format;
    wchar = ((fmt == RCFD_SHIFT_C0_N0_C1_W0CHAR) || (fmt == RCFD_SHIFT_C0_N0CHAR_C1_W0CHAR));

    n_shift = S->shifts_count;
    n_in    = Sin->shifts_count;
    n_out   = Sout->shifts_count;

    /* corona copies of the data as the shift itself will see them */
    parallel_map_phase_data(rcfd, i_phase);

    /* fetch, for every incoming shift, the source cell values from the node that owns the source
       (same route as parallel_shift_data: all nodes -> node 0 -> reorder by out2in -> destination node) */
    {
        out_vals = malloc((n_out * K + 1) * sizeof(double));

        for(i = 0; i < n_out; i++){

            c0 = Sout->c0[i];

            out_vals[i * K + 0] = _C(i_layer).data[i_phase][c0][i_dc];
            out_vals[i * K + 1] = _C(i_layer).x[c0][0];
            out_vals[i * K + 2] = _C(i_layer).x[c0][1];
            out_vals[i * K + 3] = _C(i_layer).x[c0][2];
            out_vals[i * K + 4] = _C(i_layer).vof[i_prev][c0][i_phase];
            out_vals[i * K + 5] = (c0 >= n_int) ? 1.0 : 0.0;
        }

        n_send  = n_out * K;
        n_total = PRF_GISUM1(n_send);

        recvcounts = malloc(sizeof(int) * compute_node_count);
        displs     = malloc(sizeof(int) * compute_node_count);

        PRF_GATHER_INT(recvcounts, &n_send, 1, 0);

        displs[0] = 0;

        for(i_node = 1; i_node < compute_node_count; ++i_node){

            displs[i_node] = displs[i_node-1] + recvcounts[i_node-1];
        }

        if(myid == 0){

            g_out = malloc((n_total + 1) * sizeof(double));
            g_in  = malloc((n_total + 1) * sizeof(double));
        }

        PRF_GATHERV_REAL(g_out, out_vals, recvcounts, n_send, displs, 0);

        free(recvcounts);
        free(displs);

        if(myid == 0){

            n_tr = n_total / K;

            for(i = 0; i < n_tr; i++){

                memcpy(&g_in[out2in[i] * K], &g_out[i * K], K * sizeof(double));
            }

            i_global = to_node_count[0] * K;

            for(i_node = 1; i_node < compute_node_count; i_node++){

                PRF_CSEND_INT(i_node, &to_node_count[i_node], 1, node_zero);

                if(to_node_count[i_node] > 0){

                    PRF_CSEND_REAL(i_node, &g_in[i_global], to_node_count[i_node] * K, node_zero);

                    i_global += to_node_count[i_node] * K;
                }
            }

            in_vals = g_in;
            n_mine  = to_node_count[0];
        }

        if(myid > 0){

            PRF_CRECV_INT(node_zero, &n_mine, 1, node_zero);

            if(n_mine > 0){

                in_vals = malloc(n_mine * K * sizeof(double));

                PRF_CRECV_REAL(node_zero, in_vals, n_mine * K, node_zero);
            }
        }

        if(n_mine > n_in) n_mine = n_in;
    }

#define RCFD_C2C_W(L, i_) (wchar ? ((double)(L)->w0char[i_] / 255.0) : (L)->w0[i_])

    /* weight sum of every destination cell */
    {
        wsum = calloc(n_int + 1, sizeof(double));

        for(i = 0; i < n_shift; i++){ w = RCFD_C2C_W(S, i);   if(w > 0.0) wsum[S->c1[i]]   += w; }
        for(i = 0; i < n_mine;  i++){ w = RCFD_C2C_W(Sin, i); if(w > 0.0) wsum[Sin->c1[i]] += w; }
    }

    if(rcfd_c2c_pred_size < n_int){

        free(rcfd_c2c_pred);
        rcfd_c2c_pred = malloc((n_int + 1) * sizeof(double));
        rcfd_c2c_pred_size = n_int;
    }
    pred = rcfd_c2c_pred;
    for(i = 0; i < n_int; i++) pred[i] = 0.0;

    for(j = 0; j < RCFD_C2C_DIAG_NSTAT; j++) st[j] = 0.0;

    write_file = ((step >= RCFD_C2C_DIAG_FILE_FIRST) && (step <= RCFD_C2C_DIAG_FILE_LAST) && RCFD_C2C_DIAG_FILE_NODE(myid));

    if(write_file){

        sprintf(name, "rCFD_c2c_step%d_node%d.txt", step, myid);
        fo = fopen(name, "w");

        if(fo){

            fprintf(fo, "# step %d frame %d (previous frame %d) node %d phase %d, concentration = data %d\n", step, i_frame, i_prev, myid, i_phase, i_dc);
            fprintf(fo, "# class: 0 source is a cell of this node, 1 corona cell of this node, 2 cell of another node, 3 corona cell of another node\n");
            fprintf(fo, "# c1 x1 y1 z1 vof1 vof1_prev conc1_before class source_node c0 x0 y0 z0 vof0_prev conc0 w wsum\n");
        }
    }

    /* cells: mass, composition before the shift, cells not hit */
    for(i = 0; i < n_int; i++){

        m = _C(i_layer).volume[i] * _C(i_layer).vof[i_frame][i][i_phase] * rho;

        st[0] += m;
        st[1] += m * _C(i_layer).data[i_phase][i][i_dc];

        if(!(wsum[i] > 0.0)){

            st[2] += m;
            st[3] += m * _C(i_layer).data[i_phase][i][i_dc];
        }
    }

    st[4] = (double)n_shift;
    st[5] = (double)n_mine;

    /* shifts */
    for(i = 0; i < n_shift + n_mine; i++){

        if(i < n_shift){

            w = RCFD_C2C_W(S, i);
            if(!(w > 0.0)) continue;

            c0 = S->c0[i];
            c1 = S->c1[i];
            k  = (c0 < n_int) ? 0 : 1;

            c_src        = _C(i_layer).data[i_phase][c0][i_dc];
            vof_prev_src = _C(i_layer).vof[i_prev][c0][i_phase];
            for(i_dim = 0; i_dim < 3; i_dim++) x0[i_dim] = _C(i_layer).x[c0][i_dim];
            i_node = myid;
            j = k;
        }
        else{

            int ii = i - n_shift;

            w = RCFD_C2C_W(Sin, ii);
            if(!(w > 0.0)) continue;

            c0 = Sin->c0[ii];
            c1 = Sin->c1[ii];
            k  = 2;

            c_src        = in_vals[ii * K + 0];
            x0[0]        = in_vals[ii * K + 1];
            x0[1]        = in_vals[ii * K + 2];
            x0[2]        = in_vals[ii * K + 3];
            vof_prev_src = in_vals[ii * K + 4];
            i_node = shifts_get_node0(rcfd, Sin, ii);
            j = (in_vals[ii * K + 5] > 0.5) ? 3 : 2;
        }

        m = _C(i_layer).volume[c1] * _C(i_layer).vof[i_frame][c1][i_phase] * rho;
        f = m * w / wsum[c1];

        d = 0.0;
        for(i_dim = 0; i_dim < 3; i_dim++){ dx = x0[i_dim] - _C(i_layer).x[c1][i_dim]; d += dx * dx; }
        d = sqrt(d);

        for(b = 0; b < RCFD_C2C_DIAG_NBIN - 1; b++){ if(d < bin_edge[b]) break; }

        {
            double *q = &st[6 + k * RCFD_C2C_DIAG_NCL];

            q[b]                          += f;
            q[RCFD_C2C_DIAG_NBIN + b]     += f * c_src;
            q[2 * RCFD_C2C_DIAG_NBIN]     += f * (x0[1] - _C(i_layer).x[c1][1]);
            q[2 * RCFD_C2C_DIAG_NBIN + 1] += 1.0;
        }

        if(d > md[k]) md[k] = d;

        pred[c1] += w / wsum[c1] * c_src;

        if(fo){

            fprintf(fo, "%d %.5f %.5f %.5f %.4e %.4e %.5f %d %d %d %.5f %.5f %.5f %.4e %.5f %.5f %.5f\n",
                    c1, _C(i_layer).x[c1][0], _C(i_layer).x[c1][1], _C(i_layer).x[c1][2],
                    _C(i_layer).vof[i_frame][c1][i_phase], _C(i_layer).vof[i_prev][c1][i_phase],
                    _C(i_layer).data[i_phase][c1][i_dc],
                    j, i_node, c0, x0[0], x0[1], x0[2], vof_prev_src, c_src, w, wsum[c1]);
        }
    }

    if(fo) fclose(fo);

    /* cells not hit keep their value (or are filled by shift_fill): no prediction for them */
    for(i = 0; i < n_int; i++){ if(!(wsum[i] > 0.0)) pred[i] = -1.0; }

    rcfd_c2c_pred_valid = 1;
    rcfd_c2c_pred_data  = i_dc;

    /* print, node by node */
    for(i_node = 0; i_node < compute_node_count; i_node++){

        double a_cl, b_cl, after;

        for(j = 0; j < RCFD_C2C_DIAG_NSTAT; j++) g[j] = PRF_GRSUM1((myid == i_node) ? st[j] : 0.0);
        for(k = 0; k < 3; k++)                   gmd[k] = PRF_GRHIGH1((myid == i_node) ? md[k] : 0.0);

        if(myid != 0) continue;
        if(!(g[0] > 0.0)) continue;

        after = g[3];
        for(k = 0; k < 3; k++) for(b = 0; b < RCFD_C2C_DIAG_NBIN; b++) after += g[6 + k * RCFD_C2C_DIAG_NCL + RCFD_C2C_DIAG_NBIN + b];

        report_print0(rcfd, "[c2c] step %d frame %d node %2d: solids %.4f kg, concentration %.4f before and %.4f after the shift (expected from the shift lists); %.3f %% of the solids not hit by any shift (concentration %.3f)",
                      step, i_frame, i_node, g[0], g[1] / g[0], after / g[0], 100.0 * g[2] / g[0], (g[2] > 0.0) ? g[3] / g[2] : 0.0);

        for(k = 0; k < 3; k++){

            double *q = &g[6 + k * RCFD_C2C_DIAG_NCL];

            a_cl = 0.0; b_cl = 0.0;
            for(b = 0; b < RCFD_C2C_DIAG_NBIN; b++){ a_cl += q[b]; b_cl += q[RCFD_C2C_DIAG_NBIN + b]; }

            len = snprintf(line, sizeof(line), "[c2c]   node %2d from %s: %8.0f shifts, %5.1f %% of the weight, concentration %.3f, mean dy %+6.2f mm, largest distance %7.1f mm | weight by distance 0-5 / 5-10 / 10-20 / 20-50 / >50 mm:",
                           i_node, class_name[k], q[2 * RCFD_C2C_DIAG_NBIN + 1], 100.0 * a_cl / g[0], (a_cl > 0.0) ? b_cl / a_cl : 0.0,
                           (a_cl > 0.0) ? 1000.0 * q[2 * RCFD_C2C_DIAG_NBIN] / a_cl : 0.0, 1000.0 * gmd[k]);

            for(b = 0; b < RCFD_C2C_DIAG_NBIN; b++) len += snprintf(line + len, sizeof(line) - len, " %5.1f", 100.0 * q[b] / g[0]);

            len += snprintf(line + len, sizeof(line) - len, " %% | concentration:");

            for(b = 0; b < RCFD_C2C_DIAG_NBIN; b++){

                if(q[b] > 0.0) len += snprintf(line + len, sizeof(line) - len, " %.3f", q[RCFD_C2C_DIAG_NBIN + b] / q[b]);
                else           len += snprintf(line + len, sizeof(line) - len, "   -  ");
            }

            report_print0(rcfd, "%s", line);
        }
    }

#undef RCFD_C2C_W

    free(out_vals);
    free(wsum);

    if(myid == 0){ free(g_out); free(g_in); }
    else if(in_vals){ free(in_vals); }
#endif
}

/* after shift_data and shift_fill: does the result agree with what the lists said? */
static void rcfd_c2c_diag_check(rCFD_type *rcfd, int i_phase)
{
#if RP_NODE
    int     i, n_int, n_bad = 0, i_layer = 0;
    double  d, d_max = 0.0;

    if(!rcfd_c2c_pred_valid) return;

    rcfd_c2c_pred_valid = 0;

    n_int = _CELL_DICT(i_layer).number_of_int_cells;

    if(n_int > rcfd_c2c_pred_size) return;

    for(i = 0; i < n_int; i++){

        if(rcfd_c2c_pred[i] < 0.0) continue;

        d = fabs(_C(i_layer).data[i_phase][i][rcfd_c2c_pred_data] - rcfd_c2c_pred[i]);

        if(d > d_max) d_max = d;
        if(d > 1.0e-6) n_bad++;
    }

    d_max = PRF_GRHIGH1(d_max);
    n_bad = PRF_GISUM1(n_bad);

    report_print0(rcfd, "[c2c] step %d: shift result against the value expected from the shift lists: largest difference %.3e, cells with a difference above 1e-6: %d",
                  rcfd->Solver->global_run_counter, d_max, n_bad);
#endif
}

/* ------------------------------------------------------------------------------
 * TUHH 11/2026: stage budget of a run time step (diagnostics, no effect on data).
 *
 * For a phase with two or more concentration data the species mass held by the
 * internal cells of every compute node is taken before and after each stage of
 * the time step (cell shifts, face swaps, drift, balance correction,
 * concentration-sum control). Printed per stage and data:
 *   global  : change of the global species mass in percent of the inventory
 *             (a stage that conserves the species prints 0),
 *   per node: change of the species mass on each node in 0.1 % of the inventory
 *             (what the stage moved into or out of each partition),
 * and the range of the concentration sum after the stage.
 * Reported for the first RCFD_STAGE_DIAG_FIRST steps and then every
 * RCFD_STAGE_DIAG_EVERY-th step; set both to 0 to switch off.
 * ------------------------------------------------------------------------------ */
#ifndef RCFD_STAGE_DIAG_FIRST
#define RCFD_STAGE_DIAG_FIRST 12
#endif
#ifndef RCFD_STAGE_DIAG_EVERY
#define RCFD_STAGE_DIAG_EVERY 40
#endif
#define RCFD_STAGE_DIAG_MAX_NODES 128
#define RCFD_STAGE_DIAG_MAX_DATA  4

static double rcfd_stage_prev[RCFD_STAGE_DIAG_MAX_DATA][RCFD_STAGE_DIAG_MAX_NODES];

static int rcfd_stage_diag_active(rCFD_type *rcfd, int i_phase)
{
#if RP_NODE
    int step = rcfd->Solver->global_run_counter;

    if(rcfd->Phase_Dict[i_phase].number_of_concentration_data < 2) return 0;
    if(compute_node_count > RCFD_STAGE_DIAG_MAX_NODES) return 0;
    if(rcfd->Solver->current_layer != 0) return 0;

    if(step < RCFD_STAGE_DIAG_FIRST) return 1;
    if((RCFD_STAGE_DIAG_EVERY > 0) && ((step % RCFD_STAGE_DIAG_EVERY) == 0)) return 1;
#endif
    return 0;
}

static void rcfd_stage_diag(rCFD_type *rcfd, int i_phase, const char *stage)
{
#if RP_NODE
    int     i_cell, i_data, i_node, n_data, n_conc;
    int     i_layer = 0, i_frame, int_cell_count;
    int     n_low = 0, n_out = 0, is_start;
    double  m[RCFD_STAGE_DIAG_MAX_DATA], cur[RCFD_STAGE_DIAG_MAX_DATA][RCFD_STAGE_DIAG_MAX_NODES];
    double  cell_mass, sum_c, c, sum_min = 1.0e300, sum_max = -1.0e300, inventory, total;
    char    line[2048];
    int     len;

    if(!rcfd_stage_diag_active(rcfd, i_phase)) return;

    i_frame        = rcfd->Solver->current_frame;
    int_cell_count = _CELL_DICT(i_layer).number_of_int_cells;
    n_data         = rcfd->Phase_Dict[i_phase].number_of_data;
    is_start       = (strcmp(stage, "start") == 0);

    if(n_data > RCFD_STAGE_DIAG_MAX_DATA) n_data = RCFD_STAGE_DIAG_MAX_DATA;

    for(i_data = 0; i_data < n_data; i_data++) m[i_data] = 0.0;

    for(i_cell = 0; i_cell < int_cell_count; i_cell++){

        cell_mass = _C(i_layer).volume[i_cell] * _C(i_layer).vof[i_frame][i_cell][i_phase] * rcfd->Phase_Dict[i_phase].density;

        sum_c = 0.0; n_conc = 0;

        for(i_data = 0; i_data < n_data; i_data++){

            if(rcfd->Data_Dict[i_phase][i_data].type != RCFD_DATA_CONCENTRATION) continue;

            c = _C(i_layer).data[i_phase][i_cell][i_data];

            m[i_data] += c * cell_mass;
            sum_c     += c;

            if((c < -1.0e-9) || (c > 1.0 + 1.0e-9)) n_conc++;
        }

        if(_C(i_layer).vof[i_frame][i_cell][i_phase] > rcfd->Solver_Dict->vof_cutoff_limit){

            if(sum_c < sum_min) sum_min = sum_c;
            if(sum_c > sum_max) sum_max = sum_c;
            if(sum_c < 0.5)     n_low++;
            if(n_conc > 0)      n_out++;
        }
    }

    for(i_data = 0; i_data < n_data; i_data++){

        for(i_node = 0; i_node < compute_node_count; i_node++){

            cur[i_data][i_node] = PRF_GRSUM1((myid == i_node) ? m[i_data] : 0.0);
        }
    }

    sum_min = PRF_GRLOW1(sum_min);
    sum_max = PRF_GRHIGH1(sum_max);
    n_low   = PRF_GISUM1(n_low);
    n_out   = PRF_GISUM1(n_out);

    if((myid == 0) && (!is_start)){

        for(i_data = 0; i_data < n_data; i_data++){

            if(rcfd->Data_Dict[i_phase][i_data].type != RCFD_DATA_CONCENTRATION) continue;

            inventory = rcfd->Balance[i_phase][i_data].mass_integral_target_global;
            total     = 0.0;

            for(i_node = 0; i_node < compute_node_count; i_node++) total += cur[i_data][i_node] - rcfd_stage_prev[i_data][i_node];

            if(!(inventory > 0.0)) inventory = 1.0;

            len = snprintf(line, sizeof(line), "[stage] step %d phase %d %-7s data %d: global %+8.4f %% | per node [0.1 %%]:",
                           rcfd->Solver->global_run_counter, i_phase, stage, i_data, 100.0 * total / inventory);

            for(i_node = 0; (i_node < compute_node_count) && (len < (int)sizeof(line) - 16); i_node++){

                len += snprintf(line + len, sizeof(line) - len, " %+6.1f",
                                1000.0 * (cur[i_data][i_node] - rcfd_stage_prev[i_data][i_node]) / inventory);
            }

            report_print0(rcfd, "%s", line);
        }

        report_print0(rcfd, "[stage] step %d phase %d %-7s sum of concentrations in [%.4f, %.4f]; cells with sum < 0.5: %d; cells with a concentration outside [0,1]: %d",
                      rcfd->Solver->global_run_counter, i_phase, stage, sum_min, sum_max, n_low, n_out);
    }

    for(i_data = 0; i_data < n_data; i_data++){

        for(i_node = 0; i_node < compute_node_count; i_node++){

            rcfd_stage_prev[i_data][i_node] = cur[i_data][i_node];
        }
    }
#endif
}

/**
 * @brief run rCFD timestep cycles
 *
 * @param handle    pointer to previously created rCFD instance
 */
void rcfd_run(void *handle)
{
#if RP_NODE
    int i_run, i_phase;

    rCFD_type *rcfd = (rCFD_type*)handle;

    report_start_timer(rcfd, RCFD_TIMER_RUN);

    for(i_run = 0; i_run < rcfd->Solver_Dict->number_of_runs; i_run++){

        /* Update state, layer and frame */
        {
            run_update_state(rcfd);

#if MULTI_LAYER_IMPLEMENTATION_COMPLETE
            run_update_layer(rcfd);
#endif

            report_start_timer(rcfd, RCFD_TIMER_RUN_UPDATE);
            run_update_frame(rcfd);
            post_run_update_frame(rcfd);

            run_update_balance(rcfd);
        }

        /* User access before phase loop */
        {
            report_start_timer(rcfd, RCFD_TIMER_USER);
            rcfd->rCFD_UDF->_user_access_data_before_phase_loop(handle);
            report_stop_timer(rcfd, RCFD_TIMER_USER);
        }

        for(i_phase = 0; i_phase < rcfd->Solver_Dict->number_of_phases; i_phase++){

            /* User access before shift */
            {
                report_start_timer(rcfd, RCFD_TIMER_USER);
                rcfd->rCFD_UDF->_user_access_data_before_shift(handle, i_phase);
                report_stop_timer(rcfd, RCFD_TIMER_USER);
            }

            /* Cell shifts */
            {

                rcfd_stage_diag(rcfd, i_phase, "start");

                if(rcfd->Solver_Dict->data_convection_on){

                    rcfd_c2c_diag(rcfd, i_phase);       /* TUHH 11/2026: diagnostics only */

                    report_start_timer(rcfd, RCFD_TIMER_RUN_SHIFT);

                    shift_data(rcfd, i_phase);

                    shift_fill(rcfd, i_phase);

                    post_run_shift(rcfd, i_phase);

                    rcfd_c2c_diag_check(rcfd, i_phase); /* TUHH 11/2026: diagnostics only */
                }

                rcfd_stage_diag(rcfd, i_phase, "shift");
            }

            /* Face swaps */
            {
                report_start_timer(rcfd, RCFD_TIMER_RUN_SWAP);

                if(rcfd->Solver_Dict->face_swap_diffusion_on && rcfd->Solver->current_layer == 0){

                    swap_data(rcfd, i_phase);
                }

                rcfd_stage_diag(rcfd, i_phase, "swap");

                if(rcfd->Solver_Dict->data_drifting_on && rcfd->Solver->current_layer == 0){

                    swap_drift_data(rcfd, i_phase);
                }

                rcfd_stage_diag(rcfd, i_phase, "drift");

                if(rcfd->Solver_Dict->data_binarization_on && rcfd->Solver->current_layer == 0){

                    swap_binarize_data(rcfd, i_phase);
                }

                post_run_swap(rcfd, i_phase);
            }

            /* User access after swap */
            {
                report_start_timer(rcfd, RCFD_TIMER_USER);
                rcfd->rCFD_UDF->_user_access_data_after_swap(handle, i_phase);
                report_stop_timer(rcfd, RCFD_TIMER_USER);
            }

            /* Balance correction */
            {
                report_start_timer(rcfd, RCFD_TIMER_RUN_BALANCE);

                balance_update(rcfd, i_phase);

                balance_correct(rcfd, i_phase);

                rcfd_stage_diag(rcfd, i_phase, "balance");

                balance_control_conc_sum(rcfd, i_phase);

                rcfd_stage_diag(rcfd, i_phase, "concsum");

                post_run_balance(rcfd, i_phase);
            }

        }

        /* User access after phase loop */
        {
            report_start_timer(rcfd, RCFD_TIMER_USER);
            rcfd->rCFD_UDF->_user_post(handle);
            report_stop_timer(rcfd, RCFD_TIMER_USER);
        }

        dump_data(rcfd);

        rcfd->Solver->global_run_counter++;

        rcfd->Solver->global_time += rcfd->Solver->timestep_width[rcfd->Solver->current_layer];
    }

    post_run(rcfd);
#endif
}

/**
 * @brief delete rCFD instance created by rcfd_open() or its variants
 *
 * @param handle    pointer to previously created rCFD instance
 */
void rcfd_close(void *handle)
{
    if(!handle)
        return;

    rCFD_type *rcfd = (rCFD_type*)handle;

    if(rcfd->initialized){

        int mode =  rcfd->Solver_Dict->mode;

        if(mode == RCFD_MODE_PREPARATION){

            prep_write_metadata(rcfd);
        }

        free_all(rcfd);

        post_free(rcfd, mode);
    }

#if defined(RCFD_LIB_MPI)
    if(rcfd->Parallel->finalize){

        int flag = 0;
        MPI_Initialized(&flag);

        if(flag){

            MPI_Finalized(&flag);
            if(!flag){

                free_report(rcfd);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Finalize();
            }
        }
    }
#endif

    free_rcfd(rcfd);
}

/**
 * @brief set displacement and diffusion fields
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the layer
 * @param i_phase   specifies the phase
 * @param cfd_displ (flattened) displacement field (x, y, z values for each cell)
 * @param cfd_displ_diff  (flattened) diffusion field (xx, yy, zz, xy, yz, zx values for each cell)
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_set_displacement_fields(void *handle, int i_layer, int i_phase, double *cfd_displ, double *cfd_displ_diff)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                         "[%d] rcfd_set_displacement_fields:"
                         " setting displacement fields only possible in preparation mode", myid);
            return -1;
        }

        if(    (i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            || (i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)){

            report_error(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                         "[%d] rcfd_set_displacement_fields: invalid args", myid);
            return -1;
        }

        if(cfd_displ && _C(i_layer).cfd_displ){

            memcpy(&_C(i_layer).cfd_displ[i_phase][0][0], cfd_displ, sizeof(double) * _CELL_DICT(i_layer).number_of_cells * 3);
        }

        if(cfd_displ_diff && _C(i_layer).cfd_displ_diff){

            memcpy(&_C(i_layer).cfd_displ_diff[i_phase][0][0], cfd_displ_diff, sizeof(double) * _CELL_DICT(i_layer).number_of_cells * 6);
        }

        return 0;
    }
#endif // RP_NODE

    return -1;
}

/**
 * @brief generate primary cell-to-cell shifts from continuous displacement fields,
 *        additional shifts may be generated after this function has succeeded by
 *        calling rcfd_interpolate_shifts()
 *
 * @param handle    pointer to previously created rCFD instance
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_create_shifts_from_displacement_fields(void *handle)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                         "[%d] rcfd_create_shifts_from_displacement_fields:"
                         " generation of shifts only possible in preparation mode", myid);
            return -1;
        }

        int i_layer = rcfd->Solver->current_layer;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers){

            report_error(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                         "[%d] rcfd_create_shifts_from_displacement_fields: invalid layer %d", myid, i_layer);
            return -1;
        }

        // generate primary shifts for all phases
        if(rcfd->Solver_Dict->shifts_from_displacement_fields_on){

            int i_phase;
            int phase_count = rcfd->Solver_Dict->number_of_phases;
            int primary_shift_count, global_primary_shift_count;

            for(i_phase = 0; i_phase < phase_count; i_phase++){

                primary_shift_count = prep_shift_from_displacement_fields(rcfd, i_layer, i_phase);
                global_primary_shift_count = PRF_GISUM1(primary_shift_count);

                report_info0(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                             "rcfd_create_shifts_from_displacement_fields (layer %d, phase %d): number of primary shifts: %d",
                             i_layer, i_phase, global_primary_shift_count);
            }
        }
        else{

            report_error0(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                          "rcfd_create_shifts_from_displacement_fields:"
                          " shifts_from_displacement_fields_on is turned off");
            return -1;
        }

        return 0;
    }
#endif // RP_NODE

    return -1;
}

/**
 * @brief create cell-to-cell shifts for grid interfaces (in order to prevent information locks)
 *
 * @param handle                    pointer to previously created rCFD instance
 * @param interface_face_centers    (flattened) array of (process-local) coordinates for interface face centers
 * @param interface_face_count      (process-local) number of interface face centers
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_create_shifts_from_interfaces(void   *handle,
                                       double *interface_face_centers,
                                       int     interface_face_count)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                         "[%d] rcfd_create_shifts_from_interfaces:"
                         " creation of shifts only possible in preparation mode", myid);
            return -1;
        }

        if(!rcfd->Solver_Dict->grid_interface_shifts_on){

            report_error0(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                          "rcfd_create_shifts_from_interfaces:"
                          " grid_interface_shifts_on is turned off");
            return -1;
        }

        int interface_shift_count = prep_shift_from_interfaces(rcfd,
                                                               interface_face_centers,
                                                               interface_face_count);

        int error_flag = PRF_GILOW1(interface_shift_count);

        if(error_flag < 0){

            report_error0(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                         "rcfd_create_shifts_from_interfaces: failed to create interface shifts");
        }
        else{

            int global_interface_shift_count = PRF_GISUM1(interface_shift_count);

            if(global_interface_shift_count > 0){

                report_info0(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                             "rcfd_create_shifts_from_interfaces: created %d interface shifts",
                             global_interface_shift_count);
            }
        }

        return 0;
    }
#endif // RP_NODE

    return -1;
}

/**
 * @brief Generate additional cell-to-cell shifts via interpolation for more
 *        accurate shift patterns; only available in preparation mode
 *
 * @param handle    pointer to previously created rCFD instance
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_interpolate_shifts(void *handle)
{
#if RP_NODE

    /* known limitations:
        ->  single layer
    */

    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                         "[%d] rcfd_interpolate_shifts: interpolation only possible in preparation mode", myid);
            return -1;
        }

        return prep_interpolate_shifts(rcfd);
    }
#endif

    return -1;
}

/**
 * @brief sets vof field variable
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   current grid layer
 * @param cell_vof  two-dim array of cell vof values [i_cell][i_phase]
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_set_vof(void *handle, int i_layer, double** cell_vof)
{
#if RP_NODE

    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_CREATE_SHIFT,
                         "[%d] rcfd_set_vof: only available in preparation mode", myid);
            return -1;
        }

        int i_cell, i_phase, i_frame = 0;

        int int_cell_count = _CELL_DICT(i_layer).number_of_int_cells;
        int cell_count =     _CELL_DICT(i_layer).number_of_cells;
        int phase_count =    rcfd->Solver_Dict->number_of_phases;

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            for(i_phase = 0; i_phase < phase_count; i_phase++){

                _C(i_layer).vof[i_frame][i_cell][i_phase] = cell_vof[i_cell][i_phase];
            }
        }

        double *tmp = malloc(cell_count * sizeof(double));
        memset(tmp, 0, cell_count * sizeof(double));

        for(i_phase = 0; i_phase < phase_count; i_phase++){

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                tmp[i_cell] = cell_vof[i_cell][i_phase];
            }

            parallel_communicate(rcfd, tmp, i_layer, RCFD_PARALLEL_MAP_HOST_FORMAT);

            for(i_cell = int_cell_count; i_cell < cell_count; i_cell++){

                _C(i_layer).vof[i_frame][i_cell][i_phase] = tmp[i_cell];
            }
        }

        free(tmp);

        return 0;
    }
#endif

    return -1;
}

/**
 * @brief write cell-to-cell shifts to files; only available in preparation mode
 *
 * @param handle    pointer to previously created rCFD instance
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_write_shifts(void *handle)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;
        int i_phase, i_state, i_frame, i_layer;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_SHIFT,
                         "[%d] rcfd_write_shifts: only possible in preparation mode", myid);
            return -1;
        }

        i_frame = rcfd->Solver->current_frame;

        if(i_frame >= rcfd->Solver_Dict->number_of_frames){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_SHIFT,
                         "[%d] rcfd_write_shifts: reached max. number of frames", myid);
            return -1;
        }

        {
            int stored_C2C_count;

            i_state = rcfd->Solver->current_state;
            i_layer = rcfd->Solver->current_layer;

            for(i_phase = 0; i_phase < rcfd->Solver_Dict->number_of_phases; i_phase++){

                stored_C2C_count = shift_write(rcfd, i_state, i_frame, i_phase, i_layer);
                stored_C2C_count = PRF_GISUM1(stored_C2C_count);

                report_info0(rcfd, RCFD_REPORT_CATEGORY_WRITE_SHIFT,
                             "rcfd_write_shifts: total number of shifts %d (i_state %d i_frame %d, i_phase %d i_layer %d)",
                             stored_C2C_count, i_state, i_frame, i_phase, i_layer);
            }

            rcfd->C2C_Interface->valid_shifts_count = 0;

            rcfd->Solver->current_frame++;
        }

        return 0;
    }
#endif // RP_NODE

    return -1;
}

/**
 * @brief write the recurrence matrix, recurrence path and recurrence path
 *        differences to files; only available in preparation mode
 * @param handle    pointer to previously created rCFD instance
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_write_recurrence(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_FIELDS,
                         "[%d] rcfd_write_recurrence: only possible in preparation mode", myid);
            return -1;
        }

        return prep_write_recurrence(rcfd);
    }

    return -1;
}

/**
 * @brief Write cell, face and node data to (process-specific) topology files
 * @param handle    pointer to previously created rCFD instance
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_write_topo(void *handle)
{
#if RP_NODE

    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(!rcfd->initialized){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_TOPOLOGY,
                         "[%d] rcfd_write_topo: rcfd is not initialized", myid);
            return -1;
        }

        return topo_write(rcfd);
    }
#endif

    return -1;
}

/**
 * @brief set the norms, vof, local diffusion and user-defined recurrence fields for the current frame;
 *          only available in preparation mode
 * @param handle    pointer to previously created rCFD instance
 * @param norm      array of cell norms; may be NULL if set directly or via user defined function (i.e. norms format is different from 'standard')
 * @param vof       vof field for each phase ([i_cell][i_phase]); may be NULL in case of single-phase or if set directly
 * @param local_diffusion local diffusion field ([i_phase][i_cell]); may be NULL if local diffusion is turned off or if set directly
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_set_fields(void *handle, double *norm, double **vof, double **local_diffusion)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_FIELDS,
                         "[%d] rcfd_set_fields: only possible in preparation mode", myid);
            return -1;
        }

        return prep_set_fields(rcfd, norm, vof, local_diffusion);
    }
#endif

    return -1;
}

/**
 * @brief write a database entry for norms, vof, local diffusion and user-defined recurrence fields of the current frame;
 *          only available in preparation mode
 * @param handle    pointer to previously created rCFD instance
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_write_fields(void *handle)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_FIELDS,
                         "[%d] rcfd_write_fields: only possible in preparation mode", myid);
            return -1;
        }

        return prep_write_fields(rcfd);
    }
#endif

    return -1;
}

/**
 * @brief register user defined functions
 *
 * @param handle    pointer to previously created rCFD instance
 * @param udf       name of user defined function to register to
 * @param callback  user function to be registered
 * @retval -1 Error (invalid udf string)
 * @retval 0 OK
 */
int rcfd_register_user_function(void *handle, const char *udf, void (*callback)(void))
{
    rCFD_type *rcfd = (rCFD_type*)handle;

    if(strcmp(udf, "user_set_solver_dict") == 0){

        rcfd->rCFD_UDF->_user_set_solver_dict   = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_file_dict") == 0){

        rcfd->rCFD_UDF->_user_set_file_dict     = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_phase_dict") == 0){

        rcfd->rCFD_UDF->_user_set_phase_dict    = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_prep_dict") == 0){

        rcfd->rCFD_UDF->_user_set_prep_dict     = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_norm_dict") == 0){

        rcfd->rCFD_UDF->_user_set_norm_dict     = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_rec_dict") == 0){

        rcfd->rCFD_UDF->_user_set_rec_dict      = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_data_dict") == 0){

        rcfd->rCFD_UDF->_user_set_data_dict     = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_balance_dict") == 0){

        rcfd->rCFD_UDF->_user_set_balance_dict  = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_topo_dict") == 0){

        rcfd->rCFD_UDF->_user_set_topo_dict     = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_cell_dict") == 0){

        rcfd->rCFD_UDF->_user_set_cell_dict     = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_set_face_dict") == 0){

        rcfd->rCFD_UDF->_user_set_face_dict     = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_set_node_dict") == 0){

        rcfd->rCFD_UDF->_user_set_node_dict     = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_set_dump_dict") == 0){

        rcfd->rCFD_UDF->_user_set_dump_dict     = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_pre_proc") == 0){

        rcfd->rCFD_UDF->_user_pre_proc          = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_init_cells") == 0){

        rcfd->rCFD_UDF->_user_init_cells        = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_init_faces") == 0){

        rcfd->rCFD_UDF->_user_init_faces        = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_init_cell_hosts") == 0){

        rcfd->rCFD_UDF->_user_init_cell_hosts   = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_init_nodes") == 0){

        rcfd->rCFD_UDF->_user_init_nodes        = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_init_data") == 0){

        rcfd->rCFD_UDF->_user_init_data         = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_set_layer") == 0){

        rcfd->rCFD_UDF->_user_set_layer         = (short(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_access_data_before_phase_loop") == 0){

        rcfd->rCFD_UDF->_user_access_data_before_phase_loop = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_access_data_before_shift") == 0){

        rcfd->rCFD_UDF->_user_access_data_before_shift  = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_access_data_after_swap") == 0){

        rcfd->rCFD_UDF->_user_access_data_after_swap    = (void(*)(void*,const short))callback;
    }
    else if(strcmp(udf, "user_post") == 0){

        rcfd->rCFD_UDF->_user_post              = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_norm") == 0){

        rcfd->rCFD_UDF->_user_set_norm          = (void(*)(void*))callback;
    }
    else if(strcmp(udf, "user_set_rec_vars") == 0){

        rcfd->rCFD_UDF->_user_set_rec_vars = (void(*)(void*, const short))callback;
    }
    else{

        return -1;
    }

    return 0;
}

/**
 * @brief Write a message to the node0 report output stream, depending on the message priority
 *
 * @param handle    pointer to previously created rCFD instance
 * @param priority  message priority determining visibility in output
 * @param format    string that contains the text to be written; it may contain
 *                  format specifiers that are replaced by the values specified
 *                  in subsequent additional arguments and formatted as requested
 */
void rcfd_report0(void* handle, int priority, const char *format, ...)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(priority < 0 || priority >= RCFD_REPORT_MAX){

            priority = RCFD_REPORT_UNKNOWN;
        }

        rcfd->Report->number_of_reports[RCFD_REPORT_CATEGORY_USER][priority]++;

        if(rcfd->Report->priority > priority)
            return;

        va_list args;
        va_start(args, format);
        report_write_line0(rcfd, priority, vform(format, args));
        va_end(args);
    }
}

/**
 * @brief set integer dict property
 *
 * @param handle    pointer to previously created rCFD instance
 * @param dict      name of the dictonary holding the property
 * @param i_layer   specifies the layer if the property is per-layer
 * @param i_phase   specifies the phase if the property is per-phase
 * @param i_data    specifies the data index if the property is per-data
 * @param property  name of the property to be set
 * @param value     integer value to assign to the property
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_set_dict_property_i(void *handle, const char *dict, int i_layer, int i_phase, int i_data, const char *property, int value)
{
    if(handle){

        if(strcmp(dict, "File_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "integer");
            return -1;
        }
        else if(strcmp(dict, "Solver_Dict") == 0){

            return set_solver_dict_property_i(handle, property, value);
        }
        else if(strcmp(dict, "Phase_Dict") == 0){

            return set_phase_dict_property_i(handle, i_phase, property, value);
        }
        else if(strcmp(dict, "Prep_Dict") == 0){

            return set_prep_dict_property_i(handle, i_phase, property, value);
        }
        else if(strcmp(dict, "Norm_Dict") == 0){

            return set_norm_dict_property_i(handle, property, value);
        }
        else if(strcmp(dict, "Rec_Dict") == 0){

            return set_rec_dict_property_i(handle, property, value);
        }
        else if(strcmp(dict, "Data_Dict") == 0){

            return set_data_dict_property_i(handle, i_phase, i_data, property, value);
        }
        else if(strcmp(dict, "Balance_Dict") == 0){

            return set_balance_dict_property_i(handle, i_phase, i_data, property, value);
        }
        else if(strcmp(dict, "Topo_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "integer");
            return -1;
        }
        else if(strcmp(dict, "Cell_Dict") == 0){

            return set_cell_dict_property_i(handle, i_layer, property, value);
        }
        else if(strcmp(dict, "Face_Dict") == 0){

            return set_face_dict_property_i(handle, i_layer, property, value);
        }
        else if(strcmp(dict, "Node_Dict") == 0){

            return set_node_dict_property_i(handle, i_layer, property, value);
        }
        else if(strcmp(dict, "Dump_Dict") == 0){

            return set_dump_dict_property_i(handle, i_phase, i_data, property, value);
        }
        else{

            report_dict_invalid((rCFD_type*)handle, dict, property);
            return -1;
        }
    }

    return -1;
}

/**
 * @brief set floating point dict property
 *
 * @param handle    pointer to previously created rCFD instance
 * @param dict      name of the dictonary holding the property
 * @param i_layer   specifies the layer if the property is per-layer
 * @param i_phase   specifies the phase if the property is per-phase
 * @param i_data    specifies the data index if the property is per-data
 * @param property  name of the property to be set
 * @param value     floating point value to assign to the property
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_set_dict_property_r(void *handle, const char *dict, int i_layer, int i_phase, int i_data, const char *property, double value)
{
    if(handle){

        if(strcmp(dict, "File_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "real");
            return -1;
        }
        else if(strcmp(dict, "Solver_Dict") == 0){

            return set_solver_dict_property_r(handle, property, value);
        }
        else if(strcmp(dict, "Phase_Dict") == 0){

            return set_phase_dict_property_r(handle, i_phase, property, value);
        }
        else if(strcmp(dict, "Prep_Dict") == 0){

            return set_prep_dict_property_r(handle, i_phase, property, value);
        }
        else if(strcmp(dict, "Norm_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "real");
            return -1;
        }
        else if(strcmp(dict, "Rec_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "real");
            return -1;
        }
        else if(strcmp(dict, "Data_Dict") == 0){

            return set_data_dict_property_r(handle, i_phase, i_data, property, value);
        }
        else if(strcmp(dict, "Balance_Dict") == 0){

            return set_balance_dict_property_r(handle, i_phase, i_data, property, value);
        }
        else if(strcmp(dict, "Topo_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "real");
            return -1;
        }
        else if(strcmp(dict, "Cell_Dict") == 0){

            return set_cell_dict_property_r(handle, i_layer, property, value);
        }
        else if(strcmp(dict, "Face_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "real");
            return -1;
        }
        else if(strcmp(dict, "Node_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "real");
            return -1;
        }
        else if(strcmp(dict, "Dump_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "real");
            return -1;
        }
        else{

            report_dict_invalid((rCFD_type*)handle, dict, property);
            return -1;
        }
    }

    return -1;
}

/**
 * @brief set 3d vector dict property
 *
 * @param handle    pointer to previously created rCFD instance
 * @param dict      name of the dictonary holding the property
 * @param i_layer   specifies the layer if the property is per-layer
 * @param i_phase   specifies the phase if the property is per-phase
 * @param i_data    specifies the data index if the property is per-data
 * @param property  name of the property to be set
 * @param value     3d array of double values to assign to the property
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_set_dict_property_r_3d(void *handle, const char *dict, int i_layer, int i_phase, int i_data, const char *property, double *value)
{
    if(handle){

        if(strcmp(dict, "Data_Dict") == 0){

            return set_data_dict_property_r_3d(handle, i_phase, i_data, property, value);
        }
        else if(strcmp(dict, "Dump_Dict") == 0){

            return set_dump_dict_property_r_3d(handle, i_layer, i_phase, i_data, property, value);
        }
        else{

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "vector");
            return -1;
        }
    }

    return -1;
}

/**
 * @brief set string dict property
 *
 * @param handle    pointer to previously created rCFD instance
 * @param dict      name of the dictonary holding the property
 * @param i_layer   specifies the layer if the property is per-layer
 * @param i_phase   specifies the phase if the property is per-phase
 * @param i_data    specifies the data index if the property is per-data
 * @param property  name of the property to be set
 * @param value     string value to assign to the property
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_set_dict_property_s(void *handle, const char *dict, int i_layer, int i_phase, int i_data, const char *property, const char *value)
{
    if(handle){

        if(strcmp(dict, "File_Dict") == 0){

            return set_file_dict_property_s(handle, property, value);
        }
        else if(strcmp(dict, "Solver_Dict") == 0){

            return set_solver_dict_property_s(handle, property, value);
        }
        else if(strcmp(dict, "Phase_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "string");
            return -1;
        }
        else if(strcmp(dict, "Prep_Dict") == 0){

            return set_prep_dict_property_s(handle, property, value);
        }
        else if(strcmp(dict, "Norm_Dict") == 0){

            return set_norm_dict_property_s(handle, property, value);
        }
        else if(strcmp(dict, "Rec_Dict") == 0){

            return set_rec_dict_property_s(handle, property, value);
        }
        else if(strcmp(dict, "Data_Dict") == 0){

            return set_data_dict_property_s(handle, i_phase, i_data, property, value);
        }
        else if(strcmp(dict, "Balance_Dict") == 0){

            return set_balance_dict_property_s(handle, i_phase, i_data, property, value);
        }
        else if(strcmp(dict, "Topo_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "string");
            return -1;
        }
        else if(strcmp(dict, "Cell_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "string");
            return -1;
        }
        else if(strcmp(dict, "Face_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "string");
            return -1;
        }
        else if(strcmp(dict, "Node_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "string");
            return -1;
        }
        else if(strcmp(dict, "Dump_Dict") == 0){

            report_dict_invalid_property_type((rCFD_type*)handle, dict, property, "string");
            return -1;
        }
        else{

            report_dict_invalid((rCFD_type*)handle, dict, property);
            return -1;
        }
    }

    return -1;
}

// getter methods
// from Solver_Dict_type
/**
 * @brief get global time step
 *
 * @param handle    pointer to previously created rCFD instance
 * @return global time step or 0 if handle is invalid
 */
double rcfd_get_global_timestep(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver_Dict->global_time_step;
    }

    return 0.0;
}

/**
 * @brief get vof cutoff limit
 *
 * @param handle    pointer to previously created rCFD instance
 * @return vof cutoff limit or 0 if handle is invalid
 */
double  rcfd_get_vof_cutoff_limit(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver_Dict->vof_cutoff_limit;
    }

    return 0.0;
}

/**
 * @brief get number of frames
 *
 * @param handle    pointer to previously created rCFD instance
 * @return number of frames or -1 if handle is invalid
 */
int rcfd_get_frame_count(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver_Dict->number_of_frames;
    }

    return -1;
}

/**
 * @brief get number of phases
 *
 * @param handle    pointer to previously created rCFD instance
 * @return number of phases or -1 if handle is invalid
 */
int rcfd_get_phase_count(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver_Dict->number_of_phases;
    }

    return -1;
}

/**
 * @brief get local diffusion format
 *
 * @param handle    pointer to previously created rCFD instance
 * @return diffusion format or -1 if handle is invalid
 */
int rcfd_get_local_diffusion_format(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver_Dict->local_diff_format;
    }

    return -1;
}

/**
 * @brief get number of timesteps per monitoring interval
 *
 * @param handle    pointer to previously created rCFD instance
 * @return monitoring interval size or -1 if handle is invalid
 */
int rcfd_get_monitoring_interval(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver_Dict->time_steps_per_monitoring_interval;
    }

    return -1;
}

/**
 * @brief check if vof is used in simulation
 *
 * @param handle    pointer to previously created rCFD instance
 * @retval 1 in case of multiple phases or single-phase vof
 * @retval 0 in case of single phase (vof is equal 1)
 */
int rcfd_is_vof_on(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return (rcfd->Solver_Dict->number_of_phases > 1 || rcfd->Solver_Dict->single_phase_vof);
    }

    return 0;
}

/**
 * @brief check if cell-to-cell shifts are created based on displacement fields
 *
 * @param handle    pointer to previously created rCFD instance
 * @return on/off status of setting or 0 if handle is invalid
 */
int rcfd_is_shifts_from_displacement_fields_on(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver_Dict->shifts_from_displacement_fields_on;
    }

    return 0;
}

/**
 * @brief check if geometry is moving
 *
 * @param handle    pointer to previously created rCFD instance
 * @return on/off status of setting or 0 if handle is invalid
 */
int rcfd_is_moving_geometry_on(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver_Dict->moving_geometry_on;
    }

    return 0;
}

/**
 * @brief check if local diffusion is turned on
 *
 * @param handle    pointer to previously created rCFD instance
 * @return on/off status of setting or 0 if handle is invalid
 */
int rcfd_is_local_diffusion_on(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver_Dict->local_diff_on;
    }

    return 0;
}

// from Cell_Dict_type
/**
 * @brief get number of (process-local) cells in specified layer
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return number of cells or -1 if handle or i_layer is invalid
 */
int rcfd_get_cell_count(void *handle, int i_layer)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return -1;

        return _CELL_DICT(i_layer).number_of_cells;
    }

    return -1;
}

/**
 * @brief get number of (process-local) interal cells in specified layer
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return number of internal cells or -1 if handle or i_layer is invalid
 */
int rcfd_get_int_cell_count(void *handle, int i_layer)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return -1;

        return _CELL_DICT(i_layer).number_of_int_cells;
    }

    return -1;
}

/**
 * @brief get number of user variables in specified layer
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return number of user variables or -1 if handle or i_layer is invalid
 */
int rcfd_get_cell_user_count(void *handle, int i_layer)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return -1;

        return _CELL_DICT(i_layer).number_of_user_vars;
    }

    return -1;
}

/**
 * @brief get number of recurrent user variables in specified layer
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return number of rec. user variables or -1 if handle or i_layer is invalid
 */
int rcfd_get_cell_rec_user_count(void *handle, int i_layer)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return -1;

        return _CELL_DICT(i_layer).number_of_rec_user_vars;
    }

    return -1;
}

// from Face_Dict_type
/**
 * @brief get number of (process-local) faces in specified layer
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return number of faces or -1 if handle or i_layer is invalid
 */
int rcfd_get_face_count(void *handle, int i_layer)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return -1;

        return _FACE_DICT(i_layer).number_of_faces;
    }

    return -1;
}

/**
 * @brief get number of (process-local) internal faces in specified layer
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return number of internal faces or -1 if handle or i_layer is invalid
 */
int rcfd_get_int_face_count(void *handle, int i_layer)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return -1;

        return _FACE_DICT(i_layer).number_of_int_faces;
    }

    return -1;
}

// from Phase_Dict_type
/**
 * @brief get number of data for the specified phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the phase
 * @return number of data or -1 if handle or i_phase is invalid
 */
int rcfd_get_phase_data_count(void *handle, int i_phase)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return -1;

        return rcfd->Phase_Dict[i_phase].number_of_data;
    }
#endif // RP_NODE

    return -1;
}

/**
 * @brief get density for the specified phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the phase
 * @return density or 0 if handle or i_phase is invalid
 */
double rcfd_get_phase_density(void *handle, int i_phase)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return 0.0;

        return rcfd->Phase_Dict[i_phase].density;
    }
#endif // RP_NODE

    return 0.0;
}

/**
 * @brief get heat capacity for the specified phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the phase
 * @return heat capacity or 0 if handle or i_phase is invalid
 */
double  rcfd_get_phase_heat_capacity(void *handle, int i_phase)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return 0.0;

        return rcfd->Phase_Dict[i_phase].heat_capacity;
    }
#endif // RP_NODE

    return 0.0;
}

/**
 * @brief get reference temperature for the specified phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the phase
 * @return reference temperature or 0 if handle or i_phase is invalid
 */
double  rcfd_get_phase_reference_temperature(void *handle, int i_phase)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return 0.0;

        return rcfd->Phase_Dict[i_phase].reference_temperature;
    }
#endif // RP_NODE

    return 0.0;
}

/**
 * @brief get number of user variables for the specified phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the phase
 * @return number of user variables or -1 if handle or i_phase is invalid or not RP_NODE
 */
int rcfd_get_phase_user_count(void *handle, int i_phase)
{
#if RP_NODE

    if(handle == NULL) return -1;

    rCFD_type *rcfd = (rCFD_type*)handle;

    if((i_phase < 0) || (i_phase >= rcfd->Solver_Dict->number_of_phases))
        return -1;

    return rcfd->Phase_Dict[i_phase].number_of_user_vars;

#endif // RP_NODE

    return -1;
}

/**
 * @brief get user vars for the specified phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the phase
 * @return pointer to array of user_vars or NULL if handle, i_phase is invalid or not RP_NODE
 */

double* rcfd_get_phase_user(void *handle, int i_phase)
{
#if RP_NODE

    if(handle == NULL) return NULL;

    rCFD_type *rcfd = (rCFD_type*)handle;

    if((i_phase < 0) || (i_phase >= rcfd->Solver_Dict->number_of_phases))
        return NULL;

    return rcfd->Phase_Dict[i_phase].user;

#endif // RP_NODE

    return NULL;
}

// from Solver_type
/**
 * @brief get current state
 *
 * @param handle    pointer to previously created rCFD instance
 * @return current state or -1 if handle is invalid
 */
int rcfd_get_current_state(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver->current_state;
    }

    return -1;
}

/**
 * @brief get upcoming state
 *
 * @param handle    pointer to previously created rCFD instance
 * @return upcoming state or -1 if handle is invalid
 */
int rcfd_get_next_state(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver->upcoming_state;
    }

    return -1;
}

/**
 * @brief get current layer
 *
 * @param handle    pointer to previously created rCFD instance
 * @return current layer or -1 if handle is invalid
 */
int rcfd_get_current_layer(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver->current_layer;
    }

    return -1;
}

/**
 * @brief get number of cell-to-cell shift frames (i.e. database entries) already written

 * @param handle    pointer to previously created rCFD instance
 * @return number of already written cell-to-cell shift frames or -1 if handle is invalid or not in preparation mode
 */
int rcfd_get_current_frame(void *handle)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver->current_frame;
    }
#endif

    return -1;
}

/**
 * @brief get global run counter
 *
 * @param handle    pointer to previously created rCFD instance
 * @return global run counter or -1 if handle is invalid
 */
int rcfd_get_run_counter(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver->global_run_counter;
    }

    return -1;
}

/**
 * @brief get global time
 *
 * @param handle    pointer to previously created rCFD instance
 * @return global time or 0 if handle is invalid
 */
double rcfd_get_time(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Solver->global_time;
    }

    return 0.0;
}

/**
 * @brief get timestep width for the specified layer
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return timestep width or 0 if handle or i_layer is invalid
 */
double rcfd_get_timestep(void *handle, int i_layer)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return 0.0;

        return rcfd->Solver->timestep_width[i_layer];
    }

    return 0.0;
}

// dealing with multiple cell layers
/**
 * @brief map cell data from first specified cell layer to second specified cell layer
 *
 * @param handle    pointer to previously created rCFD instance
 * @param from      layer to map from
 * @param to        layer to map to
 */
void rcfd_map_layer(void *handle, int from, int to)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

#if MULTI_LAYER_IMPLEMENTATION_COMPLETE
        layer_map_from_to(rcfd, from, to);
#else
        report_error0(rcfd, RCFD_REPORT_CATEGORY_NOTSET,
                      "Implementation of multiple cell layers is incomplete.");
#endif // MULTI_LAYER_IMPLEMENTATION_COMPLETE
    }
#endif // RP_NODE
}


// access methods

// indirect access to Cell_type arrays
/**
 * @brief set cell shape type and nodes of the given cell; cell shapes must be set in consecutive ascending order of cell indices;
 *          cell shapes must not be reset;
 * @param handle            pointer to previously created rCFD instance
 * @param i_layer           specifies the cell layer
 * @param i_cell            specifies the cell to define
 * @param cell_type         specifies cell shape type (one of CELL_TYPE as in VTK)
 * @param cell_node_ids     array containing the number of cell nodes followed by the node indices
 * @param poly_face_nodes   in case of polyhedral cells, array containing the number of polyfaces followed by the the node count of each face
 *                          otherwise might be NULL
 * @retval -1 Error
 * @retval 0 OK
 */
int rcfd_set_cell_shape(void *handle, int i_layer, int i_cell, int cell_type, const int *cell_node_ids, const int *poly_face_nodes)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return -1;

        if(i_cell < 0 || i_cell >= _CELL_DICT(i_layer).number_of_cells)
            return -1;


        return topo_set_cell_nodes(rcfd, i_layer, i_cell, cell_type, cell_node_ids, poly_face_nodes);
    }
#endif // RP_NODE

    return -1;
}

// direct access to MPI_Cell_type arrays
/**
 * @brief get (process-local) cell host array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return pointer to cell hosts or NULL if handle or i_layer is invalid
 */
int* rcfd_get_cell_host(void *handle, int i_layer)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        return _MPI_Cell.host_of_cell;
    }
#endif // RP_NODE

    return NULL;
}

// direct access to Cell_type arrays
/**
 * @brief get (process-local) cell center array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return pointer to cell centers or NULL if handle or i_layer is invalid
 */
double** rcfd_get_cell_center(void *handle, int i_layer)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        return _C(i_layer).x;
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) cell volume array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return pointer to cell volumes or NULL if handle or i_layer is invalid
 */
double* rcfd_get_cell_volume(void *handle, int i_layer)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        return _C(i_layer).volume;
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) cell data array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @param i_phase   specifies the desired phase
 * @return pointer to cell data or NULL if handle, i_layer or i_phase is invalid
 */
double** rcfd_get_cell_data(void *handle, int i_layer, int i_phase)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return NULL;

        return _C(i_layer).data[i_phase]; //[i_cell][i_data]
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) cell vof array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @param i_frame   specifies the rec frame
 * @return pointer to cell vof or NULL if handle, i_layer or i_frame is invalid
 */
double** rcfd_get_cell_vof(void *handle, int i_layer, int i_frame)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        if(i_frame < 0 || i_frame >= rcfd->Solver_Dict->number_of_frames)
            return NULL;

        return _C(i_layer).vof[i_frame]; //[i_cell][i_phase]
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) cell user data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return pointer to cell user data or NULL if handle or i_layer is invalid
 */
double** rcfd_get_cell_user(void *handle, int i_layer)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        return _C(i_layer).user; //[i_cell][i_user] (number_of_user_vars)
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) cell recurrent user data of a given frame
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @param i_frame   specifies the frame
 * @return pointer to cell rec. user data or NULL if handle or i_layer or i_frame is invalid
 */
double** rcfd_get_cell_rec_user(void *handle, int i_layer, int i_frame)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        if(i_frame < 0 || i_frame >= rcfd->Solver_Dict->number_of_frames)
            return NULL;

        if((rcfd->Solver_Dict->mode == RCFD_MODE_PREPARATION) && (i_frame != 0))
            return NULL;

        return _C(i_layer).rec_user[i_frame]; //[i_cell][i_rec_user] (number_of_rec_user_vars)
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get the value for a given probe (collective function, must be called from all processes)
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_probe   specifies the probe to query
 * @return probe value or zero if handle or i_probe is invalid
 */
double rcfd_get_probe_value(void *handle, int i_probe)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_probe < 0 || i_probe >= rcfd->Dump_Dict->probe_count)
            return 0.0;

        return dump_get_probe_value(rcfd, i_probe);
    }
#endif // RP_NODE

    return 0.0;
}

/**
 * @brief get distance between a cell center and given point coordinates
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @param i_cell    specifies the cell
 * @param point_x   specifies the point coordinates
 * @return distance or -1.0 if handle, i_layer or i_cell is invalid
 */
double rcfd_cell2point_distance(void *handle, int i_layer, int i_cell, double* point_x)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return -1.0;

        if(i_cell < 0 || i_cell >= _CELL_DICT(i_layer).number_of_cells)
            return -1.0;

        if(!point_x)
            return -1.0;

        return cell2point_distance(rcfd, i_layer, i_cell, point_x);
    }
#endif // RP_NODE

    return -1.0;
}

/**
 * @brief get (process-local) cfd diffusion data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @param i_phase   specifies the desired phase
 * @return pointer to cfd diffusion data or NULL if handle or i_layer or i_phase is invalid
 */
double** rcfd_get_cell_cfd_displ_diff(void *handle, int i_layer, int i_phase)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_FIELDS,
                         "[%d] rcfd_get_cell_cfd_displ_diff:"
                         " access to cfd diffusion only possible in preparation mode", myid);
            return NULL;
        }

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return NULL;

        return _C(i_layer).cfd_displ_diff[i_phase]; // [i_cell][i_dim]
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) local diffusion data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @param i_phase   specifies the desired phase
 * @return pointer to local diffusion data or NULL if handle or i_layer or i_phase is invalid
 */
double* rcfd_get_cell_local_diffusion(void *handle, int i_layer, int i_phase)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_FIELDS,
                         "[%d] rcfd_get_cell_local_diffusion:"
                         " access to local diffusion only possible in preparation mode", myid);
            return NULL;
        }

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return NULL;

        return _C(i_layer).local_diff[0][i_phase]; // [i_cell]
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) cell crossing time array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @param i_phase   specifies the desired phase
 * @return pointer to cell crossing times or NULL if handle, i_layer or i_phase is invalid
 */
double* rcfd_get_cell_crossing_time(void *handle, int i_layer, int i_phase)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return NULL;

        return _C(i_layer).crossing_time[i_phase]; // [i_cell]
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) cell average velocity array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @param i_phase   specifies the desired phase
 * @return pointer to cell average velocities or NULL if handle, i_layer or i_phase is invalid
 */
double* rcfd_get_cell_ave_velocity(void *handle, int i_layer, int i_phase)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return NULL;

        return _C(i_layer).average_velocity[i_phase]; // [i_cell]
    }
#endif // RP_NODE

    return NULL;
}

// direct access to Face_type arrays
/**
 * @brief get (process-local) face c0 array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return pointer to face c0 or NULL if handle or i_layer is invalid
 */
int* rcfd_get_face_c0(void *handle, int i_layer)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        return _F(i_layer).c0; //[i_face]
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) face c1 array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return pointer to face c1 or NULL if handle or i_layer is invalid
 */
int* rcfd_get_face_c1(void *handle, int i_layer)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        return _F(i_layer).c1; //[i_face]
    }
#endif // RP_NODE

    return NULL;
}

/**
 * @brief get (process-local) face area-vector array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return pointer to face area-vectors or NULL if handle or i_layer is invalid
 */
double** rcfd_get_face_area(void *handle, int i_layer)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        return _F(i_layer).area; //[i_face][i_dim]
    }
#endif // RP_NODE

    return NULL;
}

// direct access to Node_type arrays
/**
 * @brief get (process-local) node position array
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_layer   specifies the cell layer
 * @return pointer to node positions or NULL if handle or i_layer is invalid
 */
double** rcfd_get_node(void *handle, int i_layer)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_layer < 0 || i_layer >= rcfd->Solver_Dict->number_of_layers)
            return NULL;

        return _N(i_layer).x; //[i_node][i_dim]
    }
#endif // RP_NODE

    return NULL;
}

// direct access to Rec_type arrays
/**
 * @brief get jump state
 *
 * @param handle    pointer to previously created rCFD instance
 * @retval -1 Error
 * @retval 0 no jump
 * @retval 1 jumped
 */
int rcfd_get_rec_jumped_at_last_frame(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Rec->jumped_at_last_frame;
    }

    return -1;
}

/**
 * @brief get current sequence length
 *
 * @param handle    pointer to previously created rCFD instance
 * @return current sequence length or -1 if handle is invalid
 */
int rcfd_get_rec_sequence_length(void *handle)
{
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        return rcfd->Rec->sequence_length;
    }

    return -1;
}

// indirect access to Balance_type members
/**
 * @brief set the mass source for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @param ms        mass source value
 */
void rcfd_set_balance_mass_source_one(void *handle, int i_phase, int i_data, double ms)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return;

        rcfd->Balance[i_phase][i_data].mass_source = ms;
    }
#endif // RP_NODE
}

/**
 * @brief set the mass source for all data of the given phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param ms        mass source value
 */
void rcfd_set_balance_mass_source_all(void *handle, int i_phase, double ms)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        {
            int i_data;

            for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

                rcfd->Balance[i_phase][i_data].mass_source = ms;
            }
        }
    }
#endif // RP_NODE
}

/**
 * @brief set the mass source for each data of the given phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param ms        mass source values for each data
 */
void rcfd_set_balance_mass_source_each(void *handle, int i_phase, double *ms)
{
#if RP_NODE
    if(handle && ms){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        {
            int i_data;

            for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

                rcfd->Balance[i_phase][i_data].mass_source = ms[i_data];
            }
        }
    }
#endif // RP_NODE
}

/**
 * @brief modify the mass source for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @param dms       increment value to add to specified mass source
 */
void rcfd_modify_balance_mass_source_one(void *handle, int i_phase, int i_data, double dms)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return;

        rcfd->Balance[i_phase][i_data].mass_source += dms;
    }
#endif // RP_NODE
}

/**
 * @brief modify all mass sources for the given phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param dms       increment value to add to all data mass sources
 */
void rcfd_modify_balance_mass_source_all(void *handle, int i_phase, double dms)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        {
            int i_data;

            for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

                rcfd->Balance[i_phase][i_data].mass_source += dms;
            }
        }
    }
#endif // RP_NODE
}

/**
 * @brief modify each mass source for the given phase
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param dms       increment values to add to each data mass source
 */
void rcfd_modify_balance_mass_source_each(void *handle, int i_phase, double *dms)
{
#if RP_NODE
    if(handle && dms){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        {
            int i_data;

            for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

                rcfd->Balance[i_phase][i_data].mass_source += dms[i_data];
            }
        }
    }
#endif // RP_NODE
}

/**
 * @brief set the (process-local) mass integral for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @param mi        (process-local) mass integral value
 */
void rcfd_set_balance_mass_integral(void *handle, int i_phase, int i_data, double mi)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return;

        rcfd->Balance[i_phase][i_data].mass_integral = mi;
    }
#endif // RP_NODE
}

/**
 * @brief set the (process-local) mass integral target for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @param mit       (process-local) mass integral target value
 */
void rcfd_set_balance_mass_integral_target(void *handle, int i_phase, int i_data, double mit)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return;

        rcfd->Balance[i_phase][i_data].mass_integral_target = mit;
    }
#endif // RP_NODE
}

/**
 * @brief set the global mass integral for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @param mig       global mass integral value
 */
void rcfd_set_balance_mass_integral_global(void *handle, int i_phase, int i_data, double mig)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return;

        rcfd->Balance[i_phase][i_data].mass_integral_global = mig;
    }
#endif // RP_NODE
}

/**
 * @brief set the global mass integral target for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @param mitg      global mass integral target value
 */
void rcfd_set_balance_mass_integral_target_global(void *handle, int i_phase, int i_data, double mitg)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return;

        rcfd->Balance[i_phase][i_data].mass_integral_target_global = mitg;
    }
#endif // RP_NODE
}

/**
 * @brief get the (process-local) mass integral for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @return (process-local) mass integral value or 0 if handle is invalid
 */
double rcfd_get_balance_mass_integral(void *handle, int i_phase, int i_data)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return 0.0;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return 0.0;

        return rcfd->Balance[i_phase][i_data].mass_integral;
    }
#endif // RP_NODE

    return 0.0;
}

/**
 * @brief get the (process-local) mass integral target for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @return (process-local) mass integral target value or 0 if handle is invalid
 */
double rcfd_get_balance_mass_integral_target(void *handle, int i_phase, int i_data)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return 0.0;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return 0.0;

        return rcfd->Balance[i_phase][i_data].mass_integral_target;
    }
#endif // RP_NODE

    return 0.0;
}

/**
 * @brief get the global mass integral for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @return global mass integral value or 0 if handle is invalid
 */
double rcfd_get_balance_mass_integral_global(void *handle, int i_phase, int i_data)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return 0.0;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return 0.0;

        return rcfd->Balance[i_phase][i_data].mass_integral_global;
    }
#endif // RP_NODE

    return 0.0;
}

/**
 * @brief get the global mass integral target for the given phase and data
 *
 * @param handle    pointer to previously created rCFD instance
 * @param i_phase   specifies the desired phase
 * @param i_data    specifies the desired data
 * @return global mass integral target value or 0 if handle is invalid
 */
double rcfd_get_balance_mass_integral_target_global(void *handle, int i_phase, int i_data)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(i_phase < 0 || i_phase >= rcfd->Solver_Dict->number_of_phases)
            return 0.0;

        if(i_data < 0 || i_data >= rcfd->Phase_Dict[i_phase].number_of_data)
            return 0.0;

        return rcfd->Balance[i_phase][i_data].mass_integral_target_global;
    }
#endif // RP_NODE

    return 0.0;
}

/**
 * @brief get coarsening factor for norms, i.e. the reciprocal of the fraction
 *          of cells used for the calculation of the recurrence norm
 *
 * @param handle    pointer to previously created rCFD instance
 * @return norm coarsening factor or -1 if handle is invalid
 */
int rcfd_get_norm_coarsening(void *handle)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_FIELDS,
                         "[%d] rcfd_get_norm_coarsening:"
                         " access to norms only possible in preparation mode", myid);
            return -1;
        }

        return rcfd->Norm_Dict->coarse_graining;
    }
#endif

    return -1;
}

/**
 * @brief get the norm format
 * @param handle    pointer to previously created rCFD instance
 * @return the norm format or -1 if handle is invalid or not in preparation mode
 */
int rcfd_get_norm_format(void *handle)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_FIELDS,
                         "[%d] rcfd_get_norm_format:"
                         " access to norms only possible in preparation mode", myid);
            return -1;
        }

        return rcfd->Norm_Dict->format;
    }
#endif

    return -1;
}

/**
 * @brief get the (process-local) number of norms per frame; less than or equal to the
 *          (process-local) number of internal cells, depending on the norm coarsening;
 *          only available in preparation mode
 *
 * @param handle    pointer to previously created rCFD instance
 * @return the number of norms per frame or -1 if handle is invalid or not in preparation mode
 */
int rcfd_get_norm_count(void *handle)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_FIELDS,
                         "[%d] rcfd_get_norm_count:"
                         " access to norms only possible in preparation mode", myid);
            return -1;
        }

        return rcfd->Norms->number_of_norms;
    }
#endif

    return -1;
}

/**
 * @brief get (process-local) norm array, the length of the array can be obtained via
 *          rcfd_get_norm_count(); only available in preparation mode
 *
 * @param handle    pointer to previously created rCFD instance
 * @return pointer to norms or NULL if handle is invalid or not in preparation mode
 */
double *rcfd_get_norm(void *handle)
{
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(rcfd->Solver_Dict->mode != RCFD_MODE_PREPARATION){

            report_error(rcfd, RCFD_REPORT_CATEGORY_WRITE_FIELDS,
                         "[%d] rcfd_get_norm:"
                         " access to norms only possible in preparation mode", myid);
            return NULL;
        }

        return rcfd->Norms->norm;
    }
#endif

    return NULL;
}


#if defined(RCFD_LIB_TEST_MODULES)
/**
 * @brief run basic tests for the given module
 *
 * @param handle    pointer to previously created rCFD instance
 * @param module    name of module to test
 * @return 0 on success, else non-zero error code
 * @retval -1 invalid handle or module name
 */
int rcfd_test_module(void *handle, const char *module)
{
    int success = 0;
#if RP_NODE
    if(handle){

        rCFD_type *rcfd = (rCFD_type*)handle;

        if(strcmp(module, "all") == 0){

            success |= test_memory(rcfd);
            success |= test_parallel(rcfd);
            success |= test_report(rcfd);
        }
        else if(strcmp(module, "memory") == 0){

            success |= test_memory(rcfd);
        }
        else if(strcmp(module, "parallel") == 0){

            success |= test_parallel(rcfd);
        }
        else if(strcmp(module, "report") == 0){

            success |= test_report(rcfd);
        }
        else if(strcmp(module, "topo") == 0){

            success |= test_topo(rcfd);
        }
        else if(strcmp(module, "dump") == 0){

            success |= test_dump(rcfd);
        }
        else{

            return -1;
        }
    }
    else{

        return -1;
    }
#endif // RP_NODE
    return success;
}
#endif // RCFD_LIB_TEST_MODULES

