#ifndef RCFD_USER
#define RCFD_USER

#include "udf.h"
#include "rCFD_library.h"

/* rCFD run mode in OpenFOAM: a species in the gas of an aerated stirred tank.

   Case: multiphaseEuler tutorial aeratedStirredTankMRF, phases (gas liquid), gas fed through
   the patch "inlet" (sparger), leaving through "outlet".

   The gas carries one species, "tracer" (mass fraction in the gas, 0 .. 1). Gas that enters
   through the inlet between _t_feed_start and _t_feed_end has the value _inlet_value; all
   gas in the tank at the start has 0. The liquid has one species that is not used.

   After a change run   rcfd_compile run   in the case directory, then for example
       mpirun -np 8 rCFD/rcfd_run -steps 2000 -snap 20

   Output: snapshots with  c0.gas  (the species),  vof.gas  and  vf_tracer  (= species x gas
   volume fraction, zero where there is no gas), and the monitor rCFD_gas_tracer_monitor.out.
*/

#if 1

    /* names of phases: the order of   phases (gas liquid);   in system/functions */
    enum{
        gas,
        liquid,
        Phase_name_count
    };

    enum{
        tracer,
        Gas_data_name_count
    };

    enum{
        liquid_unused,
        Liquid_data_name_count
    };

/* the values of the preparation (rCFD/rCFD_user_prep.h) */
#define _global_time_step   0.025       /* frame interval [s] */
#define _number_of_frames   40

#define _density_gas        1.184
#define _density_liquid     997.0

/* feed of the species */
#define _inlet_patch        "inlet"
#define _outlet_patch       "outlet"
#define _inlet_value        1.0
#define _t_feed_start       0.0         /* [s] rCFD time */
#define _t_feed_end         1.0e30      /* a pulse: e.g. 2.0 */

#define _gas_swap_diffusion 0.0         /* face-swap diffusion of the species; 0: convection by the shifts only */

    static int Monitor_initialized = 0;

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

        /* one species per phase: nothing to normalise */
        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "control_conc_sum_on", 0);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "data_drifting_on", 0);
    }

    void user_set_phase_dict(void *handle)
    {
#if RP_NODE
        rcfd_set_dict_property_i(handle, "Phase_Dict", -1, gas   , -1, "number_of_data", Gas_data_name_count);
        rcfd_set_dict_property_i(handle, "Phase_Dict", -1, liquid, -1, "number_of_data", Liquid_data_name_count);

        rcfd_set_dict_property_r(handle, "Phase_Dict", -1, gas   , -1, "density", _density_gas);
        rcfd_set_dict_property_r(handle, "Phase_Dict", -1, liquid, -1, "density", _density_liquid);
#endif
    }

    void user_set_rec_dict(void *handle)
    {
        rcfd_set_dict_property_s(handle, "Rec_Dict", -1, -1, -1, "format", "replay");
    }

    void user_set_data_dict(void *handle)
    {
#if RP_NODE
        rcfd_set_dict_property_s(handle, "Data_Dict", -1, gas   , tracer       , "type", "concentration");
        rcfd_set_dict_property_r(handle, "Data_Dict", -1, gas   , tracer       , "face_swap_diffusion", _gas_swap_diffusion);

        rcfd_set_dict_property_s(handle, "Data_Dict", -1, liquid, liquid_unused, "type", "concentration");
        rcfd_set_dict_property_r(handle, "Data_Dict", -1, liquid, liquid_unused, "face_swap_diffusion", 0.0);
#endif
    }

    void user_set_balance_dict(void *handle)
    {
#if RP_NODE
        rcfd_set_dict_property_i(handle, "Balance_Dict", -1, gas   , tracer       , "max_correction_loops", 10);
        rcfd_set_dict_property_i(handle, "Balance_Dict", -1, liquid, liquid_unused, "max_correction_loops", 10);
#endif
    }

#endif

#if 1

    void user_init_data(void *handle, const short i_layer)
    {
#if RP_NODE
        int     i_cell;
        int     int_cell_count;
        double  **data_gas = NULL, **data_liquid = NULL;

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        data_gas    = rcfd_get_cell_data(handle, i_layer, gas   );
        data_liquid = rcfd_get_cell_data(handle, i_layer, liquid);

        if(!data_gas || !data_liquid){

            rcfd_report0(handle, RCFD_REPORT_ERROR, "user_init_data: var == NULL");

            return;
        }

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            data_gas[i_cell][tracer] = 0.0;

            data_liquid[i_cell][liquid_unused] = 0.0;
        }

        rcfd_report0(handle, RCFD_REPORT_INFO, "user_init_data");
#endif
    }

#endif

#if 1

    /* Fixed value in the cells at the inlet: the counterpart of a boundary condition.
       The cells at an inlet take their content from themselves in the shift, so the value
       set here is what the gas carries into the tank. */
    static void set_inlet_cells(void *handle)
    {
#if RP_NODE
        int         i, n;
        const int   *cells = NULL;
        double      **data_gas = NULL;
        double      current_time = rcfd_get_time(handle);

        if((current_time < _t_feed_start) || (current_time > _t_feed_end)) return;

        data_gas = rcfd_get_cell_data(handle, rcfd_get_current_layer(handle), gas);

        n = rcfd_patch_faces(handle, _inlet_patch, &cells, NULL, NULL);

        if(!data_gas || n <= 0) return;

        for(i = 0; i < n; i++){

            data_gas[cells[i]][tracer] = _inlet_value;
        }
#endif
    }

    void user_access_data_before_shift(void *handle, const short i_phase)
    {
#if RP_NODE
        if(i_phase == gas) set_inlet_cells(handle);
#endif
    }

    /* The species enters and leaves the tank, so its mass is not constant: the target of the
       mass balance is set to what is in the tank now. Without this the balance correction
       would pull the mass back to the value at the start. */
    void user_access_data_after_swap(void *handle, const short i_phase)
    {
#if RP_NODE
        int     i_layer, i_frame, i_cell;
        int     int_cell_count;
        double  mass = 0.0, mass_global;
        double  *cell_volume = NULL;
        double  **data_gas = NULL, **cell_vof = NULL;

        if(i_phase != gas) return;

        set_inlet_cells(handle);

        i_layer = rcfd_get_current_layer(handle);
        i_frame = rcfd_get_current_frame(handle);

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        cell_volume = rcfd_get_cell_volume(handle, i_layer);
        data_gas    = rcfd_get_cell_data(handle, i_layer, gas);
        cell_vof    = rcfd_get_cell_vof(handle, i_layer, i_frame);

        if(!cell_volume || !data_gas || !cell_vof) return;

        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            mass += data_gas[i_cell][tracer] * cell_volume[i_cell] * cell_vof[i_cell][gas] * _density_gas;
        }

        mass_global = PRF_GRSUM1(mass);

        rcfd_set_balance_mass_integral_target(handle, gas, tracer, mass);
        rcfd_set_balance_mass_integral_target_global(handle, gas, tracer, mass_global);
#endif
    }

#endif

#if 1

    void user_post(void *handle)
    {
#if RP_NODE
        int         i_layer, i_frame, i_cell, i, n;
        int         int_cell_count;
        const int   *cells = NULL;
        double      current_time, cell_mass_gas;
        double      mass_tracer = 0.0, mass_gas = 0.0, out_tracer = 0.0, out_gas = 0.0;
        double      *cell_volume = NULL, *vf_tracer = NULL;
        double      **data_gas = NULL, **cell_vof = NULL;

        FILE        *f_out = NULL;
        char        filename[80] = "rCFD_gas_tracer_monitor.out";

        i_layer = rcfd_get_current_layer(handle);
        i_frame = rcfd_get_current_frame(handle);

        int_cell_count = rcfd_get_int_cell_count(handle, i_layer);

        current_time = rcfd_get_time(handle);

        cell_volume = rcfd_get_cell_volume(handle, i_layer);
        data_gas    = rcfd_get_cell_data(handle, i_layer, gas);
        cell_vof    = rcfd_get_cell_vof(handle, i_layer, i_frame);

        if((data_gas == NULL)||(cell_volume == NULL)||(cell_vof == NULL)){

            rcfd_report0(handle, RCFD_REPORT_ERROR, "user_post: var == NULL");

            return;
        }

        /* species and gas in the tank */
        for(i_cell = 0; i_cell < int_cell_count; i_cell++){

            cell_mass_gas = cell_volume[i_cell] * cell_vof[i_cell][gas] * _density_gas;

            mass_gas    += cell_mass_gas;
            mass_tracer += cell_mass_gas * data_gas[i_cell][tracer];
        }

        /* species in the gas of the cells at the outlet */
        n = rcfd_patch_faces(handle, _outlet_patch, &cells, NULL, NULL);

        for(i = 0; i < n; i++){

            cell_mass_gas = cell_volume[cells[i]] * cell_vof[cells[i]][gas] * _density_gas;

            out_gas    += cell_mass_gas;
            out_tracer += cell_mass_gas * data_gas[cells[i]][tracer];
        }

        mass_gas    = PRF_GRSUM1(mass_gas);
        mass_tracer = PRF_GRSUM1(mass_tracer);
        out_gas     = PRF_GRSUM1(out_gas);
        out_tracer  = PRF_GRSUM1(out_tracer);

        if(myid == 0){

            if(Monitor_initialized == 0){

                f_out = fopen(filename, "w");

                if(f_out) fprintf(f_out, "# time [s], species mass in the tank [kg], gas mass in the tank [kg], mean species fraction in the gas, species fraction at the outlet\n");

                Monitor_initialized = 1;
            }
            else{

                f_out = fopen(filename, "a");
            }

            if(f_out){

                fprintf(f_out, "%e %e %e %e %e\n", current_time, mass_tracer, mass_gas,
                        (mass_gas > 0.0) ? mass_tracer / mass_gas : 0.0,
                        (out_gas  > 0.0) ? out_tracer  / out_gas  : 0.0);

                fclose(f_out);
            }
        }

        /* in Fluent: a user-defined memory; here: a named array, written with the snapshots.
           Species times gas volume fraction: zero where there is no gas. */
        vf_tracer = rcfd_user_field(handle, "vf_tracer");

        if(vf_tracer){

            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                vf_tracer[i_cell] = data_gas[i_cell][tracer] * cell_vof[i_cell][gas];
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
