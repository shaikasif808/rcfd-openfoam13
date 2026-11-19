#ifndef RCFD_USER
#define RCFD_USER

#include "udf.h"
#include "rCFD_library.h"

/* rCFD run mode, OpenFOAM coupling (rcfd_foam2): layering granulation in batch mode with
   three moment fields instead of size bins.

   Use: keep your present run header (cp rCFD/rCFD_user_run.h rCFD/rCFD_user_run_mixing.h), copy this
   file to   <case>/rCFD/rCFD_user_run.h,   then in the case folder
       rcfd_compile run
       mpirun -np N rCFD/rcfd_run -steps <steps> -snap <interval>
   The database is that of the mixing run (two phases); nothing has to be recorded again.

   WHAT IS CARRIED
   The particles of a cell have diameters d. With x = (d - d_ref) / d_ref the solid phase
   carries three fields per cell, each a number-weighted mean over the particles of the cell:
       X1 = <x>,  X2 = <x^2>,  X3 = <x^3>
   The shifts of rCFD mix cells linearly, and these means mix linearly as well, so the
   existing transport is used unchanged (data of type "generic").

   BALANCE
   Transport only moves particles: the three means over ALL particles change by growth alone.
   The shifts do not keep them exactly, and the balance correction of the library treats every
   field on its own, which destroys the relation between them (tested: the spread explodes).
   The library correction is therefore off for these fields and this header corrects them
   together in every step: all diameters are stretched and moved by one linear map
   x -> a x + b so that mean and variance over all particles are on target, then X3 is moved so
   that the particle volume, and with it the bed mass, is on target. The monitor reports how
   large this correction is.

   GROWTH
   Spray lands on particle surface. A cell in a spray cone receives the share
   w * (surface in the cell) of its nozzle, w = spray weight of the cell (Gaussian across the
   cone, 1/z^2 along it). Every particle of the cell then grows by the same diameter step D:
       X1 += D,  X2 += 2 D X1 + D^2,  X3 += 3 D X2 + 3 D^2 X1 + D^3      (exact, no size grid)
   D is solved from the deposited volume of the cell.

   BED MASS
   The solid field of the database only says where the particles are. The particle number is
   constant; the bed mass follows from the fields:
       M = M_core + rho_layer * N * pi/6 * ( <d^3> - <d^3>_start )

   NOT IN THIS HEADER: discharge and recycle (continuous mode), size segregation (drift),
   agglomeration, breakage, nucleation. The distribution must stay narrow and unimodal.

   Output: rCFD_granulation_monitor.out, rCFD_granulation_chambers.out and the snapshot fields
   d_mean_mm, d_sigma_um, spray_weight.
*/

/* ------------------------------------------------------------------------------------------
   settings
   ------------------------------------------------------------------------------------------ */
#if 1

    enum{ gas, solid, Phase_name_count };

    enum{ gas_mix, Gas_data_name_count };

    enum{ x1, x2, x3, Solid_data_name_count };

    /* database (the values of rCFD_user_prep.h) */
#define _global_time_step       0.005       /* frame interval [s] */
#define _number_of_frames       2000
#define _density_gas            1.2

    /* particles */
#define _density_core           1040.0      /* seed particles (gamma alumina) [kg/m3]; the solid of the database */
#define _density_layer          1440.0      /* dry layer material (sodium benzoate) [kg/m3] */
#define _layer_porosity         0.0         /* pore share of the layer; layer volume = mass / (density (1 - porosity)) */

#define _d_seed_mean            1.84e-3     /* number mean of the seeds [m] (CAMSIZER sample at t = 0) */
#define _d_seed_sigma           0.050e-3    /* standard deviation of the seeds [m] (span 0.07 at 1.84 mm) */
#define _d_ref                  1.80e-3     /* reference of x = (d - d_ref) / d_ref: the diameter of the database [m] */

    /* spray */
#define _spray_rate_kg_per_h    30.0        /* solution, all active nozzles together [kg/h] */
#define _spray_solid_content    0.30        /* mass fraction of solid in the solution */
#define _deposition_efficiency  0.97        /* share of the sprayed solid that ends on the particles */

#define _number_of_nozzles      4
    /* nozzle tips [m] in the mesh of the GF25 case (y up, z along the apparatus) */
    static const double Nozzle_tip[_number_of_nozzles][3] = {
        {0.0, 0.08, 0.125},
        {0.0, 0.08, 0.375},
        {0.0, 0.08, 0.625},
        {0.0, 0.08, 0.875}};
    static const int Nozzle_active[_number_of_nozzles] = {0, 1, 1, 0};    /* the two centre nozzles spray */
    static const double Nozzle_axis[3] = {0.0, 1.0, 0.0};                 /* upwards */

#define _spray_half_angle_deg   20.0
#define _spray_height           0.1725      /* length of the spray cone [m] */
#define _nozzle_orifice_radius  0.005
#define _atomization_pressure   2.0         /* [bar] */
#define _alpha_0                0.55        /* width of the Gaussian: sigma = alpha_0 / (1 + alpha_beta p) * cone radius */
#define _alpha_beta             0.12
#define _z_ref                  0.01        /* reference distance of the 1/z^2 decay [m] */
#define _z_min_clip             0.005
#define _spray_weight_cap       1.0
#define _spray_w_threshold      1.0e-4      /* cells below this weight are outside the cones */

    /* numerics */
#define _face_swap_diffusion    0.1         /* as in the mixing header */
#define _max_diameter_step      2.0e-6      /* limit of the diameter step of a cell per rCFD step [m]. What a cell
                                               cannot take goes to the other cells of the same cone; only if all of
                                               them are at the limit it is counted as not deposited (monitor) */
#define _growth_acceleration    1.0         /* > 1: spray rate times this factor, process time = factor x run time.
                                               The mean size is unaffected. The spread that the growth adds is too
                                               large by about this factor (variance); use 1 for width predictions. */

    /* monitors */
#define _monitor_interval       20          /* steps */
#define _number_of_chambers     4
#define _chamber_axis           2           /* z */
#define _chamber_origin         0.0
#define _chamber_length         0.25        /* [m] */

#define loop_data               for(i_data = 0; i_data < Solid_data_name_count; i_data++)
#define _pi                     3.14159265358979323846

    static double   *Spray_w = NULL;        /* spray weight per cell */
    static int      *Spray_nozzle = NULL;   /* nozzle of the cell, -1: none */

    static double   Seed_moment[3];         /* X1, X2, X3 of the seeds */
    static double   Target[3];              /* means of X1, X2, X3 over all particles as they must be */
    static double   Correction_shift = 0.0; /* sum of b of the correction x -> a x + b so far (in units of d_ref) */
    static double   Correction_last[3];     /* a - 1, b and the shift of X3 of the last step */
    static double   Deposited_mass = 0.0;   /* dry layer mass on the particles so far [kg] */
    static double   Lost_mass = 0.0;        /* sprayed solid that did not reach particles [kg] */
    static double   Limited_mass = 0.0;     /* layer mass put on particles at the limit _max_diameter_step [kg] */
    static double   Clipped_moment_sum = 0.0;
    static int      Granulation_initialized = 0;
    static long     Step_count = 0;
    static double   Solid_volume_start = 0.0;   /* solid volume of the database at the first step [m3]: fixes the
                                                   particle number and the core mass. The solid volume of the frames
                                                   varies a little; the number of particles must not. */

    static double layer_density(void)
    {
        return _density_layer * (1.0 - _layer_porosity);
    }

    /* <d^3> / d_ref^3 from the three fields */
    static double cube_mean(double X1, double X2, double X3)
    {
        return 1.0 + 3.0 * X1 + 3.0 * X2 + X3;
    }

    /* spray weight of a position and its nozzle (-1: outside all cones) */
    static double spray_weight(const double *pos, int *nozzle)
    {
        int     i, k, best = -1;
        double  w_best = 0.0;
        double  r2_cone = _nozzle_orifice_radius + _spray_height * tan(_spray_half_angle_deg * _pi / 180.0);
        double  sigma_ratio = _alpha_0 / (1.0 + _alpha_beta * _atomization_pressure);

        for(i = 0; i < _number_of_nozzles; i++){

            double d[3], z = 0.0, r = 0.0, R, sigma, z_eff, w;

            if(!Nozzle_active[i]) continue;

            for(k = 0; k < 3; k++){ d[k] = pos[k] - Nozzle_tip[i][k]; z += d[k] * Nozzle_axis[k]; }

            if(z < 0.0 || z > _spray_height) continue;

            for(k = 0; k < 3; k++){ double q = d[k] - z * Nozzle_axis[k]; r += q * q; }

            r = sqrt(r);

            R = _nozzle_orifice_radius + (z / _spray_height) * (r2_cone - _nozzle_orifice_radius);

            if(r > 1.1 * R) continue;

            sigma = sigma_ratio * R;
            if(sigma < 1.0e-6) sigma = 1.0e-6;

            z_eff = (z < _z_min_clip) ? _z_min_clip : z;

            w = (_z_ref * _z_ref) / (z_eff * z_eff) * exp(-r * r / (2.0 * sigma * sigma));

            if(w > _spray_weight_cap) w = _spray_weight_cap;

            if(w > w_best){ w_best = w; best = i; }
        }

        if(nozzle) *nozzle = (w_best > _spray_w_threshold) ? best : -1;

        return (w_best > _spray_w_threshold) ? w_best : 0.0;
    }

#endif

/* ------------------------------------------------------------------------------------------
   dictionaries
   ------------------------------------------------------------------------------------------ */
#if 1

    void user_set_solver_dict(void *handle)
    {
        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "report_priority", RCFD_REPORT_WARNING);

        rcfd_set_dict_property_r(handle, "Solver_Dict", -1, -1, -1, "global_time_step", _global_time_step);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "number_of_frames", _number_of_frames);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "number_of_phases", Phase_name_count);

        rcfd_set_dict_property_r(handle, "Solver_Dict", -1, -1, -1, "vof_cutoff_limit", 1.0e-6);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "local_diff_on", 0);

        /* the three fields are not fractions of a whole: no normalisation of their sum */
        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "control_conc_sum_on", 0);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "data_drifting_on", 0);
    }

    void user_set_phase_dict(void *handle)
    {
#if RP_NODE
        rcfd_set_dict_property_i(handle, "Phase_Dict", -1, gas  , -1, "number_of_data", Gas_data_name_count);
        rcfd_set_dict_property_i(handle, "Phase_Dict", -1, solid, -1, "number_of_data", Solid_data_name_count);

        rcfd_set_dict_property_r(handle, "Phase_Dict", -1, gas  , -1, "density", _density_gas);
        rcfd_set_dict_property_r(handle, "Phase_Dict", -1, solid, -1, "density", _density_core);
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

        rcfd_set_dict_property_s(handle, "Data_Dict", -1, gas, gas_mix, "type", "concentration");
        rcfd_set_dict_property_r(handle, "Data_Dict", -1, gas, gas_mix, "face_swap_diffusion", 1.5);

        loop_data{

            /* generic: balanced by the sum over the solid volume, i.e. over the particles */
            rcfd_set_dict_property_s(handle, "Data_Dict", -1, solid, i_data, "type", "generic");
            rcfd_set_dict_property_r(handle, "Data_Dict", -1, solid, i_data, "face_swap_diffusion", _face_swap_diffusion);
            rcfd_set_dict_property_i(handle, "Data_Dict", -1, solid, i_data, "drifting_data", 0);
        }
#endif
    }

    void user_set_balance_dict(void *handle)
    {
#if RP_NODE
        int i_data;

        /* the gas carries one constant datum that is not used */
        rcfd_set_dict_property_s(handle, "Balance_Dict", -1, gas, gas_mix, "type", "off");
        rcfd_set_dict_property_r(handle, "Balance_Dict", -1, gas, gas_mix, "accuracy_level", 1.0e30);

        /* no library correction for the three fields (see BALANCE above); with this accuracy level
           the library also does not report them */
        loop_data{

            rcfd_set_dict_property_s(handle, "Balance_Dict", -1, solid, i_data, "type", "off");
            rcfd_set_dict_property_r(handle, "Balance_Dict", -1, solid, i_data, "accuracy_level", 1.0e30);
        }
#endif
    }

#endif

/* ------------------------------------------------------------------------------------------
   start
   ------------------------------------------------------------------------------------------ */
#if 1

    void user_init_data(void *handle, const short i_layer)
    {
#if RP_NODE
        int     i_cell, int_cell_count, n_cone = 0;
        double  s, **data_solid = NULL, **data_gas = NULL, **cell_center = NULL;

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        data_gas   = rcfd_get_cell_data(handle, i_layer, gas  );
        data_solid = rcfd_get_cell_data(handle, i_layer, solid);

        cell_center = rcfd_get_cell_center(handle, i_layer);

        if(!data_gas || !data_solid || !cell_center){

            rcfd_report0(handle, RCFD_REPORT_ERROR, "user_init_data: var == NULL");

            return;
        }

        /* seeds: normal distribution around _d_seed_mean (third central moment zero) */
        s = _d_seed_sigma / _d_ref;

        Seed_moment[x1] = (_d_seed_mean - _d_ref) / _d_ref;
        Seed_moment[x2] = Seed_moment[x1] * Seed_moment[x1] + s * s;
        Seed_moment[x3] = Seed_moment[x1] * Seed_moment[x1] * Seed_moment[x1] + 3.0 * Seed_moment[x1] * s * s;

        if(Spray_w) free(Spray_w);
        if(Spray_nozzle) free(Spray_nozzle);

        Spray_w      = (double*)calloc((int_cell_count > 0) ? int_cell_count : 1, sizeof(double));
        Spray_nozzle = (int*)   calloc((int_cell_count > 0) ? int_cell_count : 1, sizeof(int));

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            data_gas[i_cell][gas_mix] = 1.0;

            data_solid[i_cell][x1] = Seed_moment[x1];
            data_solid[i_cell][x2] = Seed_moment[x2];
            data_solid[i_cell][x3] = Seed_moment[x3];

            Spray_w[i_cell] = spray_weight(cell_center[i_cell], &Spray_nozzle[i_cell]);

            if(Spray_nozzle[i_cell] >= 0) n_cone++;
        }

        n_cone = (int)(PRF_GRSUM1((double)n_cone) + 0.5);

        Target[x1] = Seed_moment[x1]; Target[x2] = Seed_moment[x2]; Target[x3] = Seed_moment[x3];

        Correction_shift = 0.0;
        Correction_last[0] = Correction_last[1] = Correction_last[2] = 0.0;

        Deposited_mass = 0.0;
        Lost_mass = 0.0;
        Limited_mass = 0.0;
        Clipped_moment_sum = 0.0;
        Step_count = 0;
        Solid_volume_start = 0.0;
        Granulation_initialized = 1;

        Message0("[user] granulation with three moment fields: seeds %.4f mm +- %.4f mm, core %.0f kg/m3, layer %.0f kg/m3\n",
                 1e3 * _d_seed_mean, 1e3 * _d_seed_sigma, (double)_density_core, layer_density());
        Message0("[user] spray %.2f kg/h solution, %.0f %% solid, efficiency %.2f: %.4f kg/h on the particles; %d cells in the spray cones\n",
                 (double)_spray_rate_kg_per_h, 100.0 * _spray_solid_content, (double)_deposition_efficiency,
                 _spray_rate_kg_per_h * _spray_solid_content * _deposition_efficiency, n_cone);

        if(_growth_acceleration != 1.0)
            Message0("[user] growth accelerated %g times: process time = %g x run time\n", (double)_growth_acceleration, (double)_growth_acceleration);

        if(n_cone == 0) Message0("[user] WARNING: no cell lies in a spray cone, check Nozzle_tip and Nozzle_active\n");

        rcfd_report0(handle, RCFD_REPORT_INFO, "user_init_data");
#endif
    }

#endif

/* ------------------------------------------------------------------------------------------
   growth
   ------------------------------------------------------------------------------------------ */
#if 1

    void user_access_data_before_phase_loop(void *handle)
    {
#if RP_NODE
        int     i_layer, i_frame, i_cell, i_data, i_nozzle, n_active = 0, int_cell_count;
        double  dt, vof_cutoff, solid_volume = 0.0, ref_volume, number_per_volume;
        double  surface_sum[_number_of_nozzles], volume_nozzle, volume_done = 0.0, volume_lost = 0.0;
        double  source[Solid_data_name_count], G[3], clipped = 0.0;
        double  scale[_number_of_nozzles], D_max, D_unit, volume_limited = 0.0;
        int     it_scale;
        double  *cell_volume = NULL, **data = NULL, **vof = NULL;

        if(!Granulation_initialized) return;

        i_layer = rcfd_get_current_layer(handle);
        i_frame = rcfd_get_current_frame(handle);

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        dt = rcfd_get_timestep(handle, i_layer);

        vof_cutoff = rcfd_get_vof_cutoff_limit(handle);

        cell_volume = rcfd_get_cell_volume(handle, i_layer);

        data = rcfd_get_cell_data(handle, i_layer, solid);

        vof = rcfd_get_cell_vof(handle, i_layer, i_frame);

        if(!cell_volume || !data || !vof) return;

        loop_data{ source[i_data] = 0.0; G[i_data] = 0.0; }

        /* means over all particles as the transport left them */
        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            double cv = cell_volume[i_cell] * vof[i_cell][solid];

            solid_volume += cv;

            loop_data{ G[i_data] += cv * data[i_cell][i_data]; }
        }

        solid_volume = PRF_GRSUM1(solid_volume);

        if(solid_volume <= 0.0) return;

        loop_data{ G[i_data] = PRF_GRSUM1(G[i_data]) / solid_volume; }

        /* correction of all three fields together: x -> a x + b puts mean and variance over all
           particles on target, then a shift of X3 puts the particle volume on target */
        {
            double a = 1.0, b, e, var_now = G[x2] - G[x1] * G[x1], var_target = Target[x2] - Target[x1] * Target[x1];

            if(var_now > 1.0e-30 && var_target > 0.0) a = sqrt(var_target / var_now);

            b = Target[x1] - a * G[x1];

            e = Target[x3] - (a * a * a * G[x3] + 3.0 * a * a * b * G[x2] + 3.0 * a * b * b * G[x1] + b * b * b);

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                double X1 = data[i_cell][x1], X2 = data[i_cell][x2], X3 = data[i_cell][x3];

                data[i_cell][x1] = a * X1 + b;
                data[i_cell][x2] = a * a * X2 + 2.0 * a * b * X1 + b * b;
                data[i_cell][x3] = a * a * a * X3 + 3.0 * a * a * b * X2 + 3.0 * a * b * b * X1 + b * b * b + e;

                /* a negative variance cannot belong to a set of particles: set it to zero and count it */
                if(data[i_cell][x2] < data[i_cell][x1] * data[i_cell][x1]){

                    clipped += (data[i_cell][x1] * data[i_cell][x1] - data[i_cell][x2]) * cell_volume[i_cell] * vof[i_cell][solid];

                    data[i_cell][x2] = data[i_cell][x1] * data[i_cell][x1];
                }
            }

            Correction_last[0] = a - 1.0; Correction_last[1] = b; Correction_last[2] = e;

            Correction_shift += b;

            Clipped_moment_sum += PRF_GRSUM1(clipped);
        }

        /* particles per solid volume of the database: the particle number is that of the seeds */
        ref_volume = _pi / 6.0 * _d_ref * _d_ref * _d_ref;

        if(Solid_volume_start <= 0.0) Solid_volume_start = solid_volume;

        /* the same number of particles in every frame: N = Solid_volume_start / (volume of a seed) */
        number_per_volume = Solid_volume_start / (solid_volume * ref_volume * cube_mean(Seed_moment[x1], Seed_moment[x2], Seed_moment[x3]));

        /* particle surface in every spray cone, weighted with the spray weight */
        for(i_nozzle = 0; i_nozzle < _number_of_nozzles; i_nozzle++){

            surface_sum[i_nozzle] = 0.0;

            if(Nozzle_active[i_nozzle]) n_active++;
        }

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            if(Spray_nozzle[i_cell] < 0 || vof[i_cell][solid] <= vof_cutoff) continue;

            /* <d^2> / d_ref^2 = 1 + 2 X1 + X2 */
            surface_sum[Spray_nozzle[i_cell]] += Spray_w[i_cell] * cell_volume[i_cell] * vof[i_cell][solid]
                                               * (1.0 + 2.0 * data[i_cell][x1] + data[i_cell][x2]);
        }

        for(i_nozzle = 0; i_nozzle < _number_of_nozzles; i_nozzle++){

            surface_sum[i_nozzle] = PRF_GRSUM1(surface_sum[i_nozzle]);
        }

        /* layer volume of one nozzle in this step [m3] */
        volume_nozzle = (n_active > 0) ? _growth_acceleration * (_spray_rate_kg_per_h / 3600.0) * _spray_solid_content
                                         * _deposition_efficiency * dt / layer_density() / n_active : 0.0;

        for(i_nozzle = 0; i_nozzle < _number_of_nozzles; i_nozzle++){

            scale[i_nozzle] = 1.0;
        }

        /* Without a limit a cell with weight w gets the diameter step  D_unit * w / surface_sum  (to first order).
           Cells above _max_diameter_step are held at the limit; scale[] raises the step of the other cells
           of the cone so that the whole spray of the nozzle is placed. */
        D_max  = _max_diameter_step * _growth_acceleration / _d_ref;
        D_unit = volume_nozzle / (number_per_volume * ref_volume * 3.0);

        for(it_scale = 0; it_scale < 6; it_scale++){

            double v_limit[_number_of_nozzles], s_free[_number_of_nozzles], change = 0.0;

            for(i_nozzle = 0; i_nozzle < _number_of_nozzles; i_nozzle++){ v_limit[i_nozzle] = 0.0; s_free[i_nozzle] = 0.0; }

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                double cv, S1, S2;

                i_nozzle = Spray_nozzle[i_cell];

                if(i_nozzle < 0 || vof[i_cell][solid] <= vof_cutoff || surface_sum[i_nozzle] <= 0.0) continue;

                cv = cell_volume[i_cell] * vof[i_cell][solid];

                S1 = 1.0 + data[i_cell][x1];
                S2 = 1.0 + 2.0 * data[i_cell][x1] + data[i_cell][x2];

                if(D_unit * Spray_w[i_cell] * scale[i_nozzle] / surface_sum[i_nozzle] > D_max)
                    v_limit[i_nozzle] += number_per_volume * cv * ref_volume * (3.0 * D_max * S2 + 3.0 * D_max * D_max * S1 + D_max * D_max * D_max);
                else
                    s_free[i_nozzle] += Spray_w[i_cell] * cv * S2;
            }

            for(i_nozzle = 0; i_nozzle < _number_of_nozzles; i_nozzle++){

                double s_new;

                if(!Nozzle_active[i_nozzle] || surface_sum[i_nozzle] <= 0.0) continue;

                v_limit[i_nozzle] = PRF_GRSUM1(v_limit[i_nozzle]);
                s_free[i_nozzle]  = PRF_GRSUM1(s_free[i_nozzle]);

                s_new = (s_free[i_nozzle] > 0.0 && volume_nozzle > v_limit[i_nozzle])
                        ? (volume_nozzle - v_limit[i_nozzle]) * surface_sum[i_nozzle] / (volume_nozzle * s_free[i_nozzle]) : 0.0;

                if(fabs(s_new - scale[i_nozzle]) > change) change = fabs(s_new - scale[i_nozzle]);

                scale[i_nozzle] = s_new;
            }

            if(change < 1.0e-9) break;
        }

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            double X1, X2, X3, S1, S2, n_cell, v_cell, target, D, f, df, cv;
            int    it;

            i_nozzle = Spray_nozzle[i_cell];

            if(i_nozzle < 0 || vof[i_cell][solid] <= vof_cutoff || surface_sum[i_nozzle] <= 0.0) continue;

            X1 = data[i_cell][x1]; X2 = data[i_cell][x2]; X3 = data[i_cell][x3];

            S1 = 1.0 + X1;                  /* <d>   / d_ref   */
            S2 = 1.0 + 2.0 * X1 + X2;       /* <d^2> / d_ref^2 */

            cv = cell_volume[i_cell] * vof[i_cell][solid];

            n_cell = number_per_volume * cv;

            /* layer volume of this cell: share of the weighted particle surface of its cone */
            v_cell = volume_nozzle * Spray_w[i_cell] * cv * S2 / surface_sum[i_nozzle] * scale[i_nozzle];

            /* all particles of the cell grow by D (in units of d_ref):
               n pi/6 d_ref^3 (3 D S2 + 3 D^2 S1 + D^3) = v_cell */
            target = v_cell / (n_cell * ref_volume);

            D = target / (3.0 * S2);

            for(it = 0; it < 4; it++){

                f  = 3.0 * D * S2 + 3.0 * D * D * S1 + D * D * D - target;
                df = 3.0 * S2 + 6.0 * D * S1 + 3.0 * D * D;

                D -= f / df;
            }

            if(D_unit * Spray_w[i_cell] * scale[i_nozzle] / surface_sum[i_nozzle] > D_max || D > D_max){

                D = D_max;

                volume_limited += n_cell * ref_volume * (3.0 * D * S2 + 3.0 * D * D * S1 + D * D * D);
            }

            volume_done += n_cell * ref_volume * (3.0 * D * S2 + 3.0 * D * D * S1 + D * D * D);

            data[i_cell][x3] = X3 + 3.0 * D * X2 + 3.0 * D * D * X1 + D * D * D;
            data[i_cell][x2] = X2 + 2.0 * D * X1 + D * D;
            data[i_cell][x1] = X1 + D;

            source[x1] += (data[i_cell][x1] - X1) * cv;
            source[x2] += (data[i_cell][x2] - X2) * cv;
            source[x3] += (data[i_cell][x3] - X3) * cv;
        }

        /* the growth of this step moves the targets */
        loop_data{ Target[i_data] += PRF_GRSUM1(source[i_data]) / solid_volume; }

        /* what was not placed (cones without particles, all cells of a cone at the limit) is lost */
        volume_done = PRF_GRSUM1(volume_done);
        volume_lost = n_active * volume_nozzle - volume_done;

        Deposited_mass += volume_done * layer_density();
        Lost_mass      += ((volume_lost > 0.0) ? volume_lost : 0.0) * layer_density();
        Limited_mass   += PRF_GRSUM1(volume_limited) * layer_density();

        Step_count++;
#endif
    }

#endif

/* ------------------------------------------------------------------------------------------
   monitors
   ------------------------------------------------------------------------------------------ */
#if 1

    void user_post(void *handle)
    {
#if RP_NODE
        static int  files_started = 0;

        int     i_layer, i_frame, i_cell, i_chamber, k, int_cell_count;
        double  sum[4 + 4 * _number_of_chambers], negative = 0.0, var_min = 1.0e30;
        double  t_run, t_process, G1, G2, G3, mean, sigma, d30, x50_3, mass_core, mass_layer, seed_cube;
        double  *cell_volume = NULL, *f_mean = NULL, *f_sigma = NULL, *f_spray = NULL;
        double  **data = NULL, **vof = NULL, **cell_center = NULL;
        FILE    *f_out = NULL;

        if(!Granulation_initialized) return;

        i_layer = rcfd_get_current_layer(handle);
        i_frame = rcfd_get_current_frame(handle);

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        cell_volume = rcfd_get_cell_volume(handle, i_layer);
        cell_center = rcfd_get_cell_center(handle, i_layer);

        data = rcfd_get_cell_data(handle, i_layer, solid);

        vof = rcfd_get_cell_vof(handle, i_layer, i_frame);

        if(!cell_volume || !cell_center || !data || !vof) return;

        /* fields of the snapshots (as C_UDMI in Fluent) */
        f_mean  = rcfd_user_field(handle, "d_mean_mm");
        f_sigma = rcfd_user_field(handle, "d_sigma_um");
        f_spray = rcfd_user_field(handle, "spray_weight");

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            double var = data[i_cell][x2] - data[i_cell][x1] * data[i_cell][x1];

            if(f_mean)  f_mean[i_cell]  = (vof[i_cell][solid] > 1.0e-3) ? 1e3 * _d_ref * (1.0 + data[i_cell][x1]) : 0.0;
            if(f_sigma) f_sigma[i_cell] = (vof[i_cell][solid] > 1.0e-3 && var > 0.0) ? 1e6 * _d_ref * sqrt(var) : 0.0;
            if(f_spray) f_spray[i_cell] = Spray_w[i_cell];
        }

        if(Step_count % _monitor_interval != 0) return;

        for(k = 0; k < 4 + 4 * _number_of_chambers; k++) sum[k] = 0.0;

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            double cv = cell_volume[i_cell] * vof[i_cell][solid];
            double var = data[i_cell][x2] - data[i_cell][x1] * data[i_cell][x1];

            if(cv <= 0.0) continue;

            sum[0] += cv;
            sum[1] += cv * data[i_cell][x1];
            sum[2] += cv * data[i_cell][x2];
            sum[3] += cv * data[i_cell][x3];

            i_chamber = (int)floor((cell_center[i_cell][_chamber_axis] - _chamber_origin) / _chamber_length);

            if(i_chamber < 0) i_chamber = 0;
            if(i_chamber >= _number_of_chambers) i_chamber = _number_of_chambers - 1;

            sum[4 + 4 * i_chamber + 0] += cv;
            sum[4 + 4 * i_chamber + 1] += cv * data[i_cell][x1];
            sum[4 + 4 * i_chamber + 2] += cv * data[i_cell][x2];
            sum[4 + 4 * i_chamber + 3] += cv * data[i_cell][x3];

            if(vof[i_cell][solid] > 1.0e-3){

                if(var < var_min) var_min = var;
                if(var < -1.0e-12) negative += 1.0;
            }
        }

        for(k = 0; k < 4 + 4 * _number_of_chambers; k++) sum[k] = PRF_GRSUM1(sum[k]);

        negative = PRF_GRSUM1(negative);
        var_min  = PRF_GRLOW1(var_min);

        if(sum[0] <= 0.0) return;

        G1 = sum[1] / sum[0]; G2 = sum[2] / sum[0]; G3 = sum[3] / sum[0];

        mean  = _d_ref * (1.0 + G1);
        sigma = _d_ref * sqrt((G2 - G1 * G1 > 0.0) ? G2 - G1 * G1 : 0.0);
        d30   = _d_ref * cbrt(cube_mean(G1, G2, G3));

        /* median of the mass distribution of a narrow, near-normal distribution */
        x50_3 = mean + 3.0 * sigma * sigma / mean;

        seed_cube = cube_mean(Seed_moment[x1], Seed_moment[x2], Seed_moment[x3]);

        /* particle number and core mass are those of the first step, not of the current frame */
        {
            double v0 = (Solid_volume_start > 0.0) ? Solid_volume_start : sum[0];

            mass_core  = _density_core * v0;
            mass_layer = layer_density() * v0 * (cube_mean(G1, G2, G3) - seed_cube) / seed_cube;
        }

        t_run     = rcfd_get_time(handle);
        t_process = t_run * _growth_acceleration;

        if(myid == 0){

            f_out = fopen("rCFD_granulation_monitor.out", files_started ? "a" : "w");

            if(f_out){

                if(!files_started)
                    fprintf(f_out, "# process time [s], run time [s], bed mass from the fields [kg], bed mass from the spray balance [kg], "
                                   "number mean d [mm], volume mean d30 [mm], standard deviation [mm], mass median x50_3 [mm] (approx.), "
                                   "span (x90-x10)/x50 (approx.), deposited [kg], not deposited [kg], cells with negative variance, "
                                   "smallest variance of x, variance set to zero so far [m3], "
                                   "mean diameter restored by the correction so far [mm], last correction: a - 1, b, shift of X3, "
                                   "deposited at the step limit [kg], solid mass of the current frame of the database [kg]\n");

                fprintf(f_out, "%e %e %e %e %e %e %e %e %e %e %e %.0f %e %e %e %e %e %e %e %e\n", t_process, t_run, mass_core + mass_layer,
                        mass_core + Deposited_mass, 1e3 * mean, 1e3 * d30, 1e3 * sigma, 1e3 * x50_3, 2.563 * sigma / x50_3,
                        Deposited_mass, Lost_mass, negative, var_min, Clipped_moment_sum,
                        1e3 * _d_ref * Correction_shift, Correction_last[0], Correction_last[1], Correction_last[2], Limited_mass, _density_core * sum[0]);

                fclose(f_out);
            }

            f_out = fopen("rCFD_granulation_chambers.out", files_started ? "a" : "w");

            if(f_out){

                if(!files_started){

                    fprintf(f_out, "# process time [s], then per chamber along axis %d (%d chambers of %g m): "
                                   "solid mass of the database [kg], number mean d [mm], standard deviation [mm]\n",
                            _chamber_axis, _number_of_chambers, (double)_chamber_length);
                }

                fprintf(f_out, "%e", t_process);

                for(i_chamber = 0; i_chamber < _number_of_chambers; i_chamber++){

                    double v = sum[4 + 4 * i_chamber], c1 = 0.0, c2 = 0.0;

                    if(v > 0.0){ c1 = sum[4 + 4 * i_chamber + 1] / v; c2 = sum[4 + 4 * i_chamber + 2] / v; }

                    fprintf(f_out, " %e %e %e", _density_core * v, 1e3 * _d_ref * (1.0 + c1),
                            1e3 * _d_ref * sqrt((c2 - c1 * c1 > 0.0) ? c2 - c1 * c1 : 0.0));
                }

                fprintf(f_out, "\n");

                fclose(f_out);
            }

            files_started = 1;
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
        REGISTER_RCFD_UDF( user_access_data_before_phase_loop );
        REGISTER_RCFD_UDF( user_post );
    }

#endif
