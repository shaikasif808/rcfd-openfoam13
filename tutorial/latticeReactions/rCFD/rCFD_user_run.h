#ifndef RCFD_USER
#define RCFD_USER

#include "udf.h"
#include "rCFD_library.h"

/* rCFD run mode in OpenFOAM: species with reactions in a steady flow through a lattice.

   Case: incompressibleFluid, one phase (fluid), channel with a Fischer-Koch lattice,
   periodic sides, patches "inlet", "outlet" and "reactingWall".

   Species (mass fractions in the fluid; the rest to 1 is inert):
       A   fed at the inlet
       B   from A in the cells at the reacting wall (fast, catalytic)
       C   from B (first order)
       D   from C (first order)
       E   from B (second order)

   One rCFD step:  shift (convection, by the library)  ->  diffusion and reactions (here).

   Diffusion and reactions are coupled: the fast reaction at the wall keeps A low in the cells
   at the wall, and that is what draws A to the wall. So they are done together, in sub-steps
   of the rCFD step: exchange between neighbour cells, then reactions, and again. With
   _diffusion_in_header 0 the library does the diffusion once per step with its face swap and
   the reactions follow once; that is simpler and under-predicts the conversion at the wall.
   With _diffusion_in_header 2 the diffusion step is implicit and contains the wall reaction:
   for meshes with thin cells at the wall, where the explicit sub-step would be far too small.

   After a change run   rcfd_compile run   in the case directory, then for example
       mpirun -np 4 rCFD/rcfd_run -steps 200 -snap 50

   Output: snapshots with c0.fluid ... c4.fluid (A ... E), the monitor
   rCFD_lattice_monitor.out (mean values in the cells at the outlet per step) and
   rCFD_lattice_profile.out (mean values per layer of cells along z, latest state).
*/

#if 1

    /* names of phases: the order of   phases (fluid);   in system/functions */
    enum{
        fluid,
        Phase_name_count
    };

    enum{
        sA,
        sB,
        sC,
        sD,
        sE,
        Fluid_data_name_count
    };

/* the values of the preparation (rCFD/rCFD_user_prep.h) */
#define _global_time_step   0.02        /* frame interval [s] */
#define _number_of_frames   2

#define _density            992.25      /* [kg/m3] */

/* diffusion of the species: D = nu (Schmidt number 1, as in the reacting CFD case).
   The face swap of the library exchanges the fraction  swap/2  of the difference of two
   neighbour cells per step; on cubes of size h this is a diffusivity  D = swap h^2 / (2 dt). */
#define _diffusivity        6.578e-7    /* [m2/s] */
#define _cell_size          2.0e-4      /* [m] uniform mesh; only for the face swap of the library (_diffusion_in_header 0) */
#define _swap_diffusion     (2.0 * _diffusivity * _global_time_step / (_cell_size * _cell_size))

/* 2: implicit diffusion with the wall reaction (for meshes with thin cells at the wall)
   1: explicit diffusion and reactions in sub-steps
   0: face swap of the library */
#define _diffusion_in_header    1
#define _max_diffusion_number   0.25    /* explicit: sub-step so that  dt_sub * sum over the faces of a cell (D A / d) / V  stays below this (must be < 1) */

/* implicit (2): the sub-step does not depend on the cell size. Each species is solved with a
   conjugate-gradient iteration; the run prints the number of iterations. */
#define _implicit_substeps      2       /* implicit steps per rCFD step */
#define _cg_tolerance           1.0e-7  /* residual / right-hand side at which the iteration stops */
#define _cg_max_iterations      2000

/* reactions [1/s]; the second-order constant is k rho / W = 1 m3/(kmol s) * 992.25 / 18 */
#define _k_AB_wall          1.0e6       /* A -> B, only in the cells at the reacting wall */
#define _k_BC               10.0        /* B -> C */
#define _k_CD               1.0         /* C -> D */
#define _k_BE               55.125      /* B + B -> E, rate k_BE Y_B^2 */
#define _reaction_substep   1.0e-3      /* [s] largest sub-step of the reaction integration */

#define _inlet_patch        "inlet"
#define _outlet_patch       "outlet"
#define _wall_patch         "reactingWall"
#define _inlet_A            1.0
#define _flow_rate          1.0e-6      /* [m3/s] through inlet and outlet, as in 0/U of the CFD */
#define _open_balance       0           /* 1: targets of the mass balances from inflow, outflow and reactions
                                           0: no correction (target = content after the step) */

/* profile along z: volume-weighted mean values in _n_layers layers of height _layer_height from _z_min */
#define _z_min              -4.0e-3
#define _layer_height       2.0e-4
#define _n_layers           100
#define _profile_interval   10          /* steps between two profile files */

    static int      Monitor_initialized = 0;
    static int      Step_count = 0;
    static int      Wall_cell_count = -1;
    static char     *Is_wall_cell = NULL;

    /* content of the domain after the last step and mean values at the outlet [kg], [-] */
    static double   Mass_last[Fluid_data_name_count];
    static double   Outlet_last[Fluid_data_name_count];
    static int      Balance_started = 0;

#endif

#if 1

    void user_set_solver_dict(void *handle)
    {
        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "report_priority", RCFD_REPORT_WARNING);

        rcfd_set_dict_property_r(handle, "Solver_Dict", -1, -1, -1, "global_time_step", _global_time_step);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "number_of_frames", _number_of_frames);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "number_of_phases", Phase_name_count);

        rcfd_set_dict_property_r(handle, "Solver_Dict", -1, -1, -1, "vof_cutoff_limit", 1.0e-6);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "local_diff_on", 0);

        /* A ... E do not add up to one (inert rest): nothing to normalise */
        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "control_conc_sum_on", 0);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "data_drifting_on", 0);
    }

    void user_set_phase_dict(void *handle)
    {
#if RP_NODE
        rcfd_set_dict_property_i(handle, "Phase_Dict", -1, fluid, -1, "number_of_data", Fluid_data_name_count);

        rcfd_set_dict_property_r(handle, "Phase_Dict", -1, fluid, -1, "density", _density);
#endif
    }

    void user_set_rec_dict(void *handle)
    {
        rcfd_set_dict_property_s(handle, "Rec_Dict", -1, -1, -1, "format", "replay");
    }

    void user_set_data_dict(void *handle)
    {
#if RP_NODE
        int i_data;

        for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

            rcfd_set_dict_property_s(handle, "Data_Dict", -1, fluid, i_data, "type", "concentration");
            rcfd_set_dict_property_r(handle, "Data_Dict", -1, fluid, i_data, "face_swap_diffusion", _diffusion_in_header ? 0.0 : _swap_diffusion);
        }
#endif
    }

    void user_set_balance_dict(void *handle)
    {
#if RP_NODE
        int i_data;

        for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

            rcfd_set_dict_property_i(handle, "Balance_Dict", -1, fluid, i_data, "max_correction_loops", 10);
        }
#endif
    }

#endif

#if 1

    void user_init_data(void *handle, const short i_layer)
    {
#if RP_NODE
        int     i_cell, i_data;
        int     int_cell_count;
        double  **data = NULL;

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        data = rcfd_get_cell_data(handle, i_layer, fluid);

        if(!data){

            rcfd_report0(handle, RCFD_REPORT_ERROR, "user_init_data: var == NULL");

            return;
        }

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                data[i_cell][i_data] = 0.0;
            }
        }

        rcfd_report0(handle, RCFD_REPORT_INFO, "user_init_data");
#endif
    }

#endif

#if 1

    /* Fixed values in the cells at the inlet: the counterpart of a boundary condition.
       The cells at an inlet take their content from themselves in the shift. */
    static void set_inlet_cells(void *handle)
    {
#if RP_NODE
        int         i, n, i_data;
        const int   *cells = NULL;
        double      **data = NULL;

        data = rcfd_get_cell_data(handle, rcfd_get_current_layer(handle), fluid);

        n = rcfd_patch_faces(handle, _inlet_patch, &cells, NULL, NULL);

        if(!data || n <= 0) return;

        for(i = 0; i < n; i++){

            for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                data[cells[i]][i_data] = 0.0;
            }

            data[cells[i]][sA] = _inlet_A;
        }
#endif
    }

    /* cells with a face on the reacting wall (the counterpart of a face loop over a wall
       thread in Fluent) */
    static void mark_wall_cells(void *handle, int int_cell_count)
    {
#if RP_NODE
        int         i, n;
        const int   *cells = NULL;

        if(Wall_cell_count >= 0) return;

        Is_wall_cell = (char*)calloc((size_t)(int_cell_count > 0 ? int_cell_count : 1), sizeof(char));

        Wall_cell_count = 0;

        n = rcfd_patch_faces(handle, _wall_patch, &cells, NULL, NULL);

        for(i = 0; i < n; i++){

            if(!Is_wall_cell[cells[i]]){ Is_wall_cell[cells[i]] = 1; Wall_cell_count++; }
        }
#endif
    }

    /* right-hand side of B, C, D, E (A is handled before) */
    static void reaction_rates(const double *y, double *r)
    {
        double rBC = _k_BC * y[sB], rCD = _k_CD * y[sC], rBE = _k_BE * y[sB] * y[sB];

        r[sA] = 0.0;
        r[sB] = -rBC - rBE;
        r[sC] =  rBC - rCD;
        r[sD] =  rCD;
        r[sE] =  rBE;
    }

    /* reactions of one cell over the time dt: the fast wall reaction exactly, the others
       with the classical Runge-Kutta scheme in sub-steps */
    static void react_cell(double *y, int at_wall, double dt)
    {
        int     i, k, n_sub;
        double  h, a, y1[Fluid_data_name_count], k1[Fluid_data_name_count], k2[Fluid_data_name_count];
        double  k3[Fluid_data_name_count], k4[Fluid_data_name_count];

        if(at_wall){

            a = y[sA] * exp(-_k_AB_wall * dt);

            y[sB] += y[sA] - a;
            y[sA]  = a;
        }

        n_sub = (int)ceil(dt / _reaction_substep - 1.0e-9);

        if(n_sub < 1) n_sub = 1;

        h = dt / (double)n_sub;

        for(k = 0; k < n_sub; k++){

            reaction_rates(y, k1);

            for(i = 0; i < Fluid_data_name_count; i++) y1[i] = y[i] + 0.5 * h * k1[i];

            reaction_rates(y1, k2);

            for(i = 0; i < Fluid_data_name_count; i++) y1[i] = y[i] + 0.5 * h * k2[i];

            reaction_rates(y1, k3);

            for(i = 0; i < Fluid_data_name_count; i++) y1[i] = y[i] + h * k3[i];

            reaction_rates(y1, k4);

            for(i = 0; i < Fluid_data_name_count; i++) y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }

    /* Diffusion and reactions over one rCFD step, in sub-steps.

       Diffusion: for every face the volume  D A / d * dt_sub  is exchanged between its two
       cells (A face area, d distance of the cell centres), the counterpart of the
       finite-volume Laplacian. The faces come from the library (rcfd_get_face_c0/c1/area);
       the faces at a processor boundary and across a periodic patch are among them, their
       second cell is a copy that rcfd_update_corona refreshes. */
    static void diffuse_and_react(void *handle, double **data, const double *cell_volume, int int_cell_count)
    {
#if RP_NODE
        static double   *conductance = NULL, *change = NULL;
        static int      n_sub = 1;

        int             i_layer = rcfd_get_current_layer(handle);
        int             n_faces = rcfd_get_face_count(handle, i_layer);
        int             *c0 = rcfd_get_face_c0(handle, i_layer), *c1 = rcfd_get_face_c1(handle, i_layer);
        double          **area = rcfd_get_face_area(handle, i_layer);
        const double    *distance = NULL;
        int             i_face, i_cell, i_data, i_sub, i, j;
        double          dt_sub, q;

        if(conductance == NULL){

            double *rate, max_rate = 0.0;

            if(rcfd_face_distance(handle, &distance) != n_faces || !c0 || !c1 || !area){

                rcfd_report0(handle, RCFD_REPORT_ERROR, "diffuse_and_react: no face data");

                return;
            }

            conductance = (double*)malloc((size_t)(n_faces > 0 ? n_faces : 1) * sizeof(double));
            change      = (double*)malloc((size_t)(int_cell_count > 0 ? int_cell_count : 1) * Fluid_data_name_count * sizeof(double));
            rate        = (double*)calloc((size_t)(int_cell_count > 0 ? int_cell_count : 1), sizeof(double));

            for(i_face = 0; i_face < n_faces; i_face++){

                conductance[i_face] = _diffusivity
                                    * sqrt(area[i_face][0] * area[i_face][0] + area[i_face][1] * area[i_face][1] + area[i_face][2] * area[i_face][2])
                                    / distance[i_face];

                rate[c0[i_face]] += conductance[i_face];

                if(c1[i_face] < int_cell_count) rate[c1[i_face]] += conductance[i_face];
            }

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                if(rate[i_cell] / cell_volume[i_cell] > max_rate) max_rate = rate[i_cell] / cell_volume[i_cell];
            }

            free(rate);

            max_rate = PRF_GRHIGH1(max_rate);

            n_sub = (int)ceil(_global_time_step * max_rate / _max_diffusion_number - 1.0e-9);

            if(n_sub < 1) n_sub = 1;

            if(myid == 0) printf("\n[user] diffusion and reactions in %d sub-steps of %g s per rCFD step\n", n_sub, _global_time_step / n_sub);
        }

        dt_sub = _global_time_step / (double)n_sub;

        for(i_sub = 0; i_sub < n_sub; i_sub++){

            rcfd_update_corona(handle, fluid);

            for(i = 0; i < int_cell_count * Fluid_data_name_count; i++) change[i] = 0.0;

            for(i_face = 0; i_face < n_faces; i_face++){

                i = c0[i_face]; j = c1[i_face];

                for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                    q = conductance[i_face] * dt_sub * (data[j][i_data] - data[i][i_data]);

                    change[i * Fluid_data_name_count + i_data] += q;

                    if(j < int_cell_count) change[j * Fluid_data_name_count + i_data] -= q;
                }
            }

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                    data[i_cell][i_data] += change[i_cell * Fluid_data_name_count + i_data] / cell_volume[i_cell];
                }
            }

            set_inlet_cells(handle);

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                react_cell(data[i_cell], Is_wall_cell[i_cell], dt_sub);
            }
        }
#endif
    }

#if RP_NODE
    /* Implicit diffusion step of one species over dt:
           V (c - c_old) / dt = sum over the faces [ D A / d (c_neighbour - c) ] - V rate c + source
       rate [1/s] takes the species out (A in the cells at the reacting wall), source [m3/s] is
       what a reaction feeds in. Cells with fixed[i] != 0 keep the value they have in x (inlet).
       The matrix is symmetric and positive definite: conjugate gradients with the diagonal as
       preconditioner. All arrays have one value per cell; x and p also for the copies of the
       cells of the neighbour processors (rcfd_update_corona_field). The iteration stops when
       the residual, weighted with the diagonal, is _cg_tolerance times that of the right-hand
       side. Returns the number of iterations. */
    static int implicit_solve(void *handle, int n_faces, const int *c0, const int *c1, const double *conductance,
                              const double *sum_conductance, const double *cell_volume, int n_int, double dt,
                              const double *rate, const double *source, const char *fixed,
                              double *x, double *r, double *z, double *p, double *q, double *diag)
    {
        int     i, j, i_face, it;
        double  loc[2], glob[2], b2, rz, rz_new, pq, alpha, beta, g;

        loc[0] = 0.0;

        for(i = 0; i < n_int; i++){

            diag[i] = cell_volume[i] / dt + sum_conductance[i] + (rate ? rate[i] * cell_volume[i] : 0.0);

            /* right-hand side in r; x holds the old values = the start of the iteration */
            r[i] = cell_volume[i] / dt * x[i] + (source ? source[i] : 0.0);

            if(!fixed[i]) loc[0] += r[i] * r[i] / diag[i];
        }

        MPI_Allreduce(loc, glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        b2 = glob[0];

        /* r = b - M x */
        rcfd_update_corona_field(handle, x);

        for(i = 0; i < n_int; i++) r[i] -= diag[i] * x[i];

        for(i_face = 0; i_face < n_faces; i_face++){

            i = c0[i_face]; j = c1[i_face]; g = conductance[i_face];

            r[i] += g * x[j];

            if(j < n_int) r[j] += g * x[i];
        }

        loc[0] = 0.0;

        for(i = 0; i < n_int; i++){

            if(fixed[i]) r[i] = 0.0;

            z[i] = r[i] / diag[i]; p[i] = z[i]; loc[0] += r[i] * z[i];
        }

        MPI_Allreduce(loc, glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        rz = glob[0];

        if(b2 <= 0.0) b2 = rz;

        for(it = 0; it < _cg_max_iterations && rz > _cg_tolerance * _cg_tolerance * b2; it++){

            /* q = M p; p is zero in the fixed cells */
            rcfd_update_corona_field(handle, p);

            for(i = 0; i < n_int; i++) q[i] = diag[i] * p[i];

            for(i_face = 0; i_face < n_faces; i_face++){

                i = c0[i_face]; j = c1[i_face]; g = conductance[i_face];

                q[i] -= g * p[j];

                if(j < n_int) q[j] -= g * p[i];
            }

            loc[0] = 0.0;

            for(i = 0; i < n_int; i++){ if(fixed[i]) q[i] = 0.0; loc[0] += p[i] * q[i]; }

            MPI_Allreduce(loc, glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            pq = glob[0];

            if(pq <= 0.0) break;

            alpha = rz / pq;

            loc[0] = 0.0;

            for(i = 0; i < n_int; i++){

                x[i] += alpha * p[i];
                r[i] -= alpha * q[i];
                z[i]  = r[i] / diag[i];

                loc[0] += r[i] * z[i];
            }

            MPI_Allreduce(loc, glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            rz_new = glob[0];
            beta   = rz_new / rz;
            rz     = rz_new;

            for(i = 0; i < n_int; i++) p[i] = z[i] + beta * p[i];
        }

        return it;
    }
#endif

    /* Diffusion and reactions over one rCFD step with implicit diffusion (_diffusion_in_header 2).
       Per sub-step: A with its sink in the cells at the reacting wall, B with what the wall
       reaction makes of A as source, C, D, E, then the reactions in the liquid in every cell.
       The wall reaction is inside the diffusion step of A, so the sub-step is not limited by
       the size of the cells at the wall. */
    static void diffuse_and_react_implicit(void *handle, double **data, const double *cell_volume, int int_cell_count)
    {
#if RP_NODE
        static double   *conductance = NULL, *sum_conductance, *rate_A, *source;
        static char     *fixed;
        static double   *x, *r, *z, *p, *q, *diag;
        static int      n_call = 0;

        int             i_layer = rcfd_get_current_layer(handle);
        int             n_faces = rcfd_get_face_count(handle, i_layer), n_cells = rcfd_get_cell_count(handle, i_layer);
        int             *c0 = rcfd_get_face_c0(handle, i_layer), *c1 = rcfd_get_face_c1(handle, i_layer);
        double          **area = rcfd_get_face_area(handle, i_layer);
        const double    *distance = NULL;
        const int       *inlet_cells = NULL;
        int             i_face, i_cell, i_data, i_sub, i, n, it[Fluid_data_name_count];
        double          dt = _global_time_step / (double)_implicit_substeps;

        if(conductance == NULL){

            size_t n_i = (size_t)(int_cell_count > 0 ? int_cell_count : 1), n_c = (size_t)(n_cells > 0 ? n_cells : 1);

            if(rcfd_face_distance(handle, &distance) != n_faces || !c0 || !c1 || !area){

                rcfd_report0(handle, RCFD_REPORT_ERROR, "diffuse_and_react_implicit: no face data");

                return;
            }

            conductance     = (double*)malloc((size_t)(n_faces > 0 ? n_faces : 1) * sizeof(double));
            sum_conductance = (double*)calloc(n_i, sizeof(double));
            rate_A          = (double*)calloc(n_i, sizeof(double));
            fixed           = (char*)  calloc(n_i, sizeof(char));
            source          = (double*)calloc(n_i, sizeof(double));
            x               = (double*)calloc(n_c, sizeof(double));
            p               = (double*)calloc(n_c, sizeof(double));
            r               = (double*)calloc(n_i, sizeof(double));
            z               = (double*)calloc(n_i, sizeof(double));
            q               = (double*)calloc(n_i, sizeof(double));
            diag            = (double*)calloc(n_i, sizeof(double));

            for(i_face = 0; i_face < n_faces; i_face++){

                conductance[i_face] = _diffusivity
                                    * sqrt(area[i_face][0] * area[i_face][0] + area[i_face][1] * area[i_face][1] + area[i_face][2] * area[i_face][2])
                                    / distance[i_face];

                sum_conductance[c0[i_face]] += conductance[i_face];

                if(c1[i_face] < int_cell_count) sum_conductance[c1[i_face]] += conductance[i_face];
            }

            /* A is taken out in the cells at the reacting wall; the cells at the inlet keep their values */
            for(i_cell = 0; i_cell < int_cell_count; i_cell++){ if(Is_wall_cell[i_cell]) rate_A[i_cell] = _k_AB_wall; }

            n = rcfd_patch_faces(handle, _inlet_patch, &inlet_cells, NULL, NULL);

            for(i = 0; i < n; i++) fixed[inlet_cells[i]] = 1;

            if(myid == 0) printf("\n[user] implicit diffusion and wall reaction, %d sub-steps of %g s per rCFD step\n", _implicit_substeps, dt);
        }

        for(i_sub = 0; i_sub < _implicit_substeps; i_sub++){

            set_inlet_cells(handle);

            for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                for(i_cell = 0; i_cell < int_cell_count; i_cell++) x[i_cell] = data[i_cell][i_data];

                if(i_data == sA){

                    it[i_data] = implicit_solve(handle, n_faces, c0, c1, conductance, sum_conductance, cell_volume, int_cell_count, dt,
                                                rate_A, NULL, fixed, x, r, z, p, q, diag);

                    /* what the wall reaction takes out of A in this sub-step goes into B [m3/s] */
                    for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                        source[i_cell] = Is_wall_cell[i_cell] ? _k_AB_wall * cell_volume[i_cell] * x[i_cell] : 0.0;
                    }
                }
                else{

                    it[i_data] = implicit_solve(handle, n_faces, c0, c1, conductance, sum_conductance, cell_volume, int_cell_count, dt,
                                                NULL, (i_data == sB) ? source : NULL, fixed, x, r, z, p, q, diag);
                }

                for(i_cell = 0; i_cell < int_cell_count; i_cell++) data[i_cell][i_data] = x[i_cell];
            }

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                react_cell(data[i_cell], 0, dt);
            }

            set_inlet_cells(handle);
        }

        if(myid == 0 && (n_call < 3 || n_call % 50 == 0)){

            printf("[user] step %d: iterations of the diffusion step for A B C D E: %d %d %d %d %d%s\n", n_call, it[sA], it[sB], it[sC], it[sD], it[sE],
                   (it[sA] >= _cg_max_iterations || it[sB] >= _cg_max_iterations) ? "  (NOT converged: raise _cg_max_iterations)" : "");
        }

        n_call++;
#endif
    }

    void user_access_data_before_shift(void *handle, const short i_phase)
    {
#if RP_NODE
        if(i_phase == fluid) set_inlet_cells(handle);
#endif
    }

    /* mean values in the cells at the outlet (volume weighted; the flow there is uniform) */
    static void outlet_mean(void *handle, double **data, const double *cell_volume, double *mean)
    {
#if RP_NODE
        int         i, n, i_data;
        const int   *cells = NULL;
        double      volume = 0.0;

        for(i_data = 0; i_data < Fluid_data_name_count; i_data++) mean[i_data] = 0.0;

        n = rcfd_patch_faces(handle, _outlet_patch, &cells, NULL, NULL);

        for(i = 0; i < n; i++){

            volume += cell_volume[cells[i]];

            for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                mean[i_data] += cell_volume[cells[i]] * data[cells[i]][i_data];
            }
        }

        volume = PRF_GRSUM1(volume);

        for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

            mean[i_data] = PRF_GRSUM1(mean[i_data]);

            if(volume > 0.0) mean[i_data] /= volume;
        }
#endif
    }

    static void domain_mass(double **data, const double *cell_volume, int int_cell_count, double *mass, double *mass_global)
    {
#if RP_NODE
        int i_cell, i_data;

        for(i_data = 0; i_data < Fluid_data_name_count; i_data++) mass[i_data] = 0.0;

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                mass[i_data] += data[i_cell][i_data] * cell_volume[i_cell] * _density;
            }
        }

        for(i_data = 0; i_data < Fluid_data_name_count; i_data++) mass_global[i_data] = PRF_GRSUM1(mass[i_data]);
#endif
    }

    /* after convection and diffusion: inlet values, reactions, and the targets of the mass
       balances.

       The shift is not conservative: the content of the domain after a step is not exactly
       the content before plus what entered minus what left. In a closed system the library
       corrects this with the mass balance. Here the system is open and reacts, so the target
       of each species is computed in this function:
           content after the last step + (inflow - outflow) dt + change by the reactions.
       The library then corrects the fields to this target. */
    void user_access_data_after_swap(void *handle, const short i_phase)
    {
#if RP_NODE
        int     i_layer, i_cell, i_data;
        int     int_cell_count;
        double  mass[Fluid_data_name_count], mass_global[Fluid_data_name_count];
        double  before[Fluid_data_name_count], before_global[Fluid_data_name_count];
        double  target_global, inflow;
        double  *cell_volume = NULL;
        double  **data = NULL;

        if(i_phase != fluid) return;

        i_layer = rcfd_get_current_layer(handle);

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        cell_volume = rcfd_get_cell_volume(handle, i_layer);
        data        = rcfd_get_cell_data(handle, i_layer, fluid);

        if(!cell_volume || !data) return;

        mark_wall_cells(handle, int_cell_count);

        set_inlet_cells(handle);

        if(!Balance_started){

            /* start: content and outlet values of the initial state */
            for(i_data = 0; i_data < Fluid_data_name_count; i_data++){ Mass_last[i_data] = 0.0; Outlet_last[i_data] = 0.0; }

            Balance_started = 1;
        }

        domain_mass(data, cell_volume, int_cell_count, before, before_global);

        if(_diffusion_in_header == 2){

            diffuse_and_react_implicit(handle, data, cell_volume, int_cell_count);
        }
        else if(_diffusion_in_header){

            diffuse_and_react(handle, data, cell_volume, int_cell_count);
        }
        else{

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                react_cell(data[i_cell], Is_wall_cell[i_cell], _global_time_step);
            }
        }

        domain_mass(data, cell_volume, int_cell_count, mass, mass_global);

        for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

            if(_open_balance){

                inflow = (i_data == sA) ? _inlet_A : 0.0;

                target_global = Mass_last[i_data]
                              + _density * _flow_rate * _global_time_step * (inflow - Outlet_last[i_data])
                              + (mass_global[i_data] - before_global[i_data]);

                if(target_global < 0.0) target_global = 0.0;
            }
            else{

                target_global = mass_global[i_data];
            }

            /* the share of this processor */
            rcfd_set_balance_mass_integral_target(handle, fluid, i_data,
                (mass_global[i_data] > 0.0) ? mass[i_data] * target_global / mass_global[i_data] : mass[i_data]);

            rcfd_set_balance_mass_integral_target_global(handle, fluid, i_data, target_global);

            Mass_last[i_data] = target_global;
        }
#endif
    }

#endif

#if 1

    void user_post(void *handle)
    {
#if RP_NODE
        int         i_layer, i_cell, i_data, i, n, k;
        int         int_cell_count;
        const int   *cells = NULL;
        double      current_time;
        double      out[Fluid_data_name_count], out_volume = 0.0, total[Fluid_data_name_count];
        double      *cell_volume = NULL;
        double      **data = NULL, **cell_center = NULL;

        FILE        *f_out = NULL;

        i_layer = rcfd_get_current_layer(handle);

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        current_time = rcfd_get_time(handle) + _global_time_step;      /* state after this step */

        cell_volume = rcfd_get_cell_volume(handle, i_layer);
        cell_center = rcfd_get_cell_center(handle, i_layer);
        data        = rcfd_get_cell_data(handle, i_layer, fluid);

        if((data == NULL)||(cell_volume == NULL)||(cell_center == NULL)){

            rcfd_report0(handle, RCFD_REPORT_ERROR, "user_post: var == NULL");

            return;
        }

        Step_count++;

        /* for the balance targets of the next step: state after the correction of this step */
        {
            double m[Fluid_data_name_count], mg[Fluid_data_name_count];

            outlet_mean(handle, data, cell_volume, Outlet_last);

            domain_mass(data, cell_volume, int_cell_count, m, mg);

            for(i_data = 0; i_data < Fluid_data_name_count; i_data++) Mass_last[i_data] = mg[i_data];
        }

        /* mean values in the cells at the outlet, and the content of the domain */
        for(i_data = 0; i_data < Fluid_data_name_count; i_data++){ out[i_data] = 0.0; total[i_data] = 0.0; }

        n = rcfd_patch_faces(handle, _outlet_patch, &cells, NULL, NULL);

        for(i = 0; i < n; i++){

            out_volume += cell_volume[cells[i]];

            for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                out[i_data] += cell_volume[cells[i]] * data[cells[i]][i_data];
            }
        }

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                total[i_data] += cell_volume[i_cell] * data[i_cell][i_data] * _density;
            }
        }

        out_volume = PRF_GRSUM1(out_volume);

        for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

            out[i_data]   = PRF_GRSUM1(out[i_data]);
            total[i_data] = PRF_GRSUM1(total[i_data]);
        }

        if(myid == 0){

            if(Monitor_initialized == 0){

                f_out = fopen("rCFD_lattice_monitor.out", "w");

                if(f_out) fprintf(f_out, "# time [s], mean A B C D E in the cells at the outlet, mass of A B C D E in the domain [kg]\n");

                Monitor_initialized = 1;
            }
            else{

                f_out = fopen("rCFD_lattice_monitor.out", "a");
            }

            if(f_out){

                fprintf(f_out, "%e", current_time);

                for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                    fprintf(f_out, " %e", (out_volume > 0.0) ? out[i_data] / out_volume : 0.0);
                }

                for(i_data = 0; i_data < Fluid_data_name_count; i_data++) fprintf(f_out, " %e", total[i_data]);

                fprintf(f_out, "\n");

                fclose(f_out);
            }
        }

        /* profile along z */
        if(Step_count % _profile_interval == 0){

            static double sum[_n_layers][Fluid_data_name_count + 1];

            for(k = 0; k < _n_layers; k++){

                for(i_data = 0; i_data <= Fluid_data_name_count; i_data++) sum[k][i_data] = 0.0;
            }

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                k = (int)floor((cell_center[i_cell][2] - _z_min) / _layer_height);

                if(k < 0 || k >= _n_layers) continue;

                for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                    sum[k][i_data] += cell_volume[i_cell] * data[i_cell][i_data];
                }

                sum[k][Fluid_data_name_count] += cell_volume[i_cell];
            }

            for(k = 0; k < _n_layers; k++){

                for(i_data = 0; i_data <= Fluid_data_name_count; i_data++){

                    sum[k][i_data] = PRF_GRSUM1(sum[k][i_data]);
                }
            }

            if(myid == 0){

                f_out = fopen("rCFD_lattice_profile.out", "w");

                if(f_out){

                    fprintf(f_out, "# time %e s: z [m], mean A B C D E in the layer of cells\n", current_time);

                    for(k = 0; k < _n_layers; k++){

                        fprintf(f_out, "%e", _z_min + (k + 0.5) * _layer_height);

                        for(i_data = 0; i_data < Fluid_data_name_count; i_data++){

                            fprintf(f_out, " %e", (sum[k][Fluid_data_name_count] > 0.0) ? sum[k][i_data] / sum[k][Fluid_data_name_count] : 0.0);
                        }

                        fprintf(f_out, "\n");
                    }

                    fclose(f_out);
                }
            }
        }
#endif
    }

#endif

    void user_register_functions(void *handle)
    {
        REGISTER_RCFD_UDF( user_set_solver_dict );
        REGISTER_RCFD_UDF( user_set_phase_dict );
        REGISTER_RCFD_UDF( user_set_rec_dict );
        REGISTER_RCFD_UDF( user_set_data_dict );
        REGISTER_RCFD_UDF( user_set_balance_dict );

        REGISTER_RCFD_UDF( user_init_data );
        REGISTER_RCFD_UDF( user_access_data_before_shift );
        REGISTER_RCFD_UDF( user_access_data_after_swap );
        REGISTER_RCFD_UDF( user_post );
    }

#endif
