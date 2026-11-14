#ifndef RCFD_USER
#define RCFD_USER

#include "udf.h"
#include "rCFD_library.h"

/* rCFD run mode on a database written inside OpenFOAM: user header.

   Same functions as the Fluent run header. After a change run   rcfd_compile run
   in the case directory, then   mpirun -np N rCFD/rcfd_run -steps 800 -snap 10

   This example: gas with one data, solids with _number_of_particle_fractions size
   classes that drift along _vert with  u_i = -k_d (d_i - d32) |d_i - d32|.
   More or fewer classes: change _number_of_particle_fractions and the two lists below.
   Another model: change the data of the phases, user_set_data_dict, user_init_data and
   user_post, as in a Fluent run header.
*/

#if 1

    /* names of phases: the order used in the preparation */
    enum{
        gas,
        solid,
        Phase_name_count
    };

    enum{
        gas_mix,
        Gas_data_name_count
    };

#define _number_of_particle_fractions   2

    static const double Diameter[_number_of_particle_fractions]      = {0.0018, 0.0033};   /* [m] */
    static const double Mass_fraction[_number_of_particle_fractions] = {0.6, 0.4};         /* initial */

#define _vert               1           /* vertical axis: 0 x, 1 y, 2 z */

#define _diam_32            0.0022      /* [m] */

#define _k_diam             2.0e5       /* drift coefficient k_d; 0: no drift */

#define _density_gas        1.2
#define _density_solid      2500.0

#define _tracer_height      -1.0        /* > 0: mixing test, class 0 below this height [m], last class above */

#define loop_diam           for(i_diam = 0; i_diam < _number_of_particle_fractions; i_diam++)

    static int Mean_height_initialized = 0;

#endif

#if 1

    void user_set_solver_dict(void *handle)
    {
        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "report_priority", RCFD_REPORT_WARNING);

        /* the values of the preparation (rCFD_user_prep.h) */
        rcfd_set_dict_property_r(handle, "Solver_Dict", -1, -1, -1, "global_time_step", 0.025);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "number_of_frames", 80);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "number_of_phases", Phase_name_count);

        rcfd_set_dict_property_r(handle, "Solver_Dict", -1, -1, -1, "vof_cutoff_limit", 1.0e-6);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "local_diff_on", 0);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "control_conc_sum_on", 1);

        /* must be 1 here if any data drifts: the drift velocities are allocated with it */
        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "data_drifting_on", (_k_diam != 0.0) ? 1 : 0);
    }

    void user_set_phase_dict(void *handle)
    {
#if RP_NODE
        rcfd_set_dict_property_i(handle, "Phase_Dict", -1, gas  , -1, "number_of_data", Gas_data_name_count);
        rcfd_set_dict_property_i(handle, "Phase_Dict", -1, solid, -1, "number_of_data", _number_of_particle_fractions);

        rcfd_set_dict_property_r(handle, "Phase_Dict", -1, gas  , -1, "density", _density_gas);
        rcfd_set_dict_property_r(handle, "Phase_Dict", -1, solid, -1, "density", _density_solid);
#endif
    }

    void user_set_rec_dict(void *handle)
    {
        rcfd_set_dict_property_s(handle, "Rec_Dict", -1, -1, -1, "format", "replay");
    }

    void user_set_data_dict(void *handle)
    {
#if RP_NODE
        int i_diam;
        double dd;
        double drift_velocity[3];

        rcfd_set_dict_property_s(handle, "Data_Dict", -1, gas  , gas_mix  , "type", "concentration");
        rcfd_set_dict_property_r(handle, "Data_Dict", -1, gas  , gas_mix  , "face_swap_diffusion", 1.5);

        drift_velocity[0] = 0.0;
        drift_velocity[1] = 0.0;
        drift_velocity[2] = 0.0;

        loop_diam{

            rcfd_set_dict_property_s(handle, "Data_Dict", -1, solid, i_diam, "type", "concentration");
            rcfd_set_dict_property_r(handle, "Data_Dict", -1, solid, i_diam, "face_swap_diffusion", 0.1);

            rcfd_set_dict_property_i(handle, "Data_Dict", -1, solid, i_diam, "drifting_data", 1);

            rcfd_set_dict_property_r(handle, "Data_Dict", -1, solid, i_diam, "d1_diam", Diameter[i_diam]);

            if(_k_diam != 0.0){

                dd = Diameter[i_diam] - _diam_32;

                drift_velocity[_vert] = - _k_diam * dd * fabs(dd);

                rcfd_set_dict_property_r_3d(handle, "Data_Dict", -1, solid, i_diam, "drift_velocity", drift_velocity);
            }
        }
#endif
    }

    void user_set_balance_dict(void *handle)
    {
#if RP_NODE
        int i_diam;

        rcfd_set_dict_property_i(handle, "Balance_Dict", -1, gas, gas_mix, "max_correction_loops", 10);

        loop_diam{

            rcfd_set_dict_property_i(handle, "Balance_Dict", -1, solid, i_diam, "max_correction_loops", 10);
        }
#endif
    }

#endif

#if 1

    void user_init_data(void *handle, const short i_layer)
    {
#if RP_NODE
        int     i_cell, i_diam;
        int     int_cell_count;
        double  **data_solid = NULL, **data_gas = NULL, **cell_center = NULL;

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        data_gas   = rcfd_get_cell_data(handle, i_layer, gas  );
        data_solid = rcfd_get_cell_data(handle, i_layer, solid);

        cell_center = rcfd_get_cell_center(handle, i_layer);

        if(!data_gas || !data_solid || !cell_center){

            rcfd_report0(handle, RCFD_REPORT_ERROR, "user_init_data: var == NULL");

            return;
        }

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            data_gas[i_cell][gas_mix] = 0.0;

            if(_tracer_height > 0.0){

                loop_diam{ data_solid[i_cell][i_diam] = 0.0; }

                if(cell_center[i_cell][_vert] < _tracer_height) data_solid[i_cell][0] = 1.0;
                else                                            data_solid[i_cell][_number_of_particle_fractions - 1] = 1.0;
            }
            else{

                loop_diam{ data_solid[i_cell][i_diam] = Mass_fraction[i_diam]; }
            }
        }

        rcfd_report0(handle, RCFD_REPORT_INFO, "user_init_data");
#endif
    }

#endif

#if 1

    void user_post(void *handle)
    {
#if RP_NODE
        int     i_layer, i_frame, i_cell, i_diam;
        int     int_cell_count;
        double  density_solid, cell_mass_solid, frac_mass, current_time;
        double  mean_height_diam[_number_of_particle_fractions];
        double  global_mass_diam[_number_of_particle_fractions];
        double  *cell_volume = NULL;
        double  **data_solid = NULL, **cell_vof = NULL, **cell_center = NULL;

        FILE    *f_out = NULL;
        char    filename_height[80] = "rCFD_mean_height_monitor.out";

        i_layer = rcfd_get_current_layer(handle);
        i_frame = rcfd_get_current_frame(handle);

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        density_solid = rcfd_get_phase_density(handle, solid);

        current_time = rcfd_get_time(handle);

        cell_volume = rcfd_get_cell_volume(handle, i_layer);
        cell_center = rcfd_get_cell_center(handle, i_layer);

        data_solid = rcfd_get_cell_data(handle, i_layer, solid);

        cell_vof = rcfd_get_cell_vof(handle, i_layer, i_frame);

        if((data_solid == NULL)||(cell_volume == NULL)||(cell_vof == NULL)||(cell_center == NULL)){

            rcfd_report0(handle, RCFD_REPORT_ERROR, "user_post: var == NULL");

            return;
        }

        loop_diam{

            mean_height_diam[i_diam] = 0.0;
            global_mass_diam[i_diam] = 0.0;
        }

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            cell_mass_solid = cell_volume[i_cell] * cell_vof[i_cell][solid] * density_solid;

            loop_diam{

                frac_mass = cell_mass_solid * data_solid[i_cell][i_diam];

                global_mass_diam[i_diam] += frac_mass;
                mean_height_diam[i_diam] += frac_mass * cell_center[i_cell][_vert];
            }
        }

        loop_diam{

            global_mass_diam[i_diam] = PRF_GRSUM1(global_mass_diam[i_diam]);
            mean_height_diam[i_diam] = PRF_GRSUM1(mean_height_diam[i_diam]);

            if(global_mass_diam[i_diam] > 1e-12){

                mean_height_diam[i_diam] /= global_mass_diam[i_diam];
            }
            else{

                mean_height_diam[i_diam] = 0.0;
            }
        }

        if(myid == 0){

            if(Mean_height_initialized == 0){

                f_out = fopen(filename_height, "w");

                if(f_out) fprintf(f_out, "# time, mean_height_frac0 .. , mass_frac0 ..\n");

                Mean_height_initialized = 1;
            }
            else{

                f_out = fopen(filename_height, "a");
            }

            if(f_out){

                fprintf(f_out, "%e", current_time);

                loop_diam{ fprintf(f_out, " %e", mean_height_diam[i_diam]); }
                loop_diam{ fprintf(f_out, " %e", global_mass_diam[i_diam]); }

                fprintf(f_out, "\n");

                fclose(f_out);
            }
        }

        /* in Fluent: C_UDMI(c, t, UDMI_vf_diam_i) = ...; here: named arrays, written with the snapshots */
        loop_diam{

            char    name[32];
            double  *vf_diam;

            sprintf(name, "vf_diam_%d", i_diam);

            vf_diam = rcfd_user_field(handle, name);

            if(vf_diam){

                for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                    vf_diam[i_cell] = data_solid[i_cell][i_diam] * cell_vof[i_cell][solid];
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
        REGISTER_RCFD_UDF( user_post );
    }

#endif
