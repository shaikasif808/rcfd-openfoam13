#ifndef RCFD_USER
#define RCFD_USER

#include "udf.h"
#include "rCFD_library.h"

/* rCFD preparation mode inside OpenFOAM: user header.

   Same functions as the Fluent prep header. After a change run   rcfd_compile prep
   in the case directory, then start foamRun.

   What is set where:
     here                    frame interval, number of frames, number of phases, interpolation,
                             vof cutoff, densities, recurrence format
     system/functions        which OpenFOAM phases are recorded (their order = phase index here)
     (entry rCFDPrep1)       and when the recording starts (recordStart)
*/

    /* user Names */
    /*************************************************************************************/

    /* names of phases: the order of   phases (air particles);   in system/functions */
    enum{
        ph_gas,
        ph_solid,

        number_of_phase_names
    };

    /* user set Dicts */
    /*************************************************************************************/

    void user_set_solver_dict(void *handle)
    {
        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "report_priority", RCFD_REPORT_WARNING);

        rcfd_set_dict_property_r(handle, "Solver_Dict", -1, -1, -1, "global_time_step", 0.025);     /* frame interval [s] */

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "number_of_frames", 80);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "number_of_phases", number_of_phase_names);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "interpolation_on", 1);

        rcfd_set_dict_property_r(handle, "Solver_Dict", -1, -1, -1, "vof_cutoff_limit", 1.0e-6);

        rcfd_set_dict_property_i(handle, "Solver_Dict", -1, -1, -1, "control_conc_sum_on", 1);

        /* shifts_from_displacement_fields_on = 1 and time_steps_per_monitoring_interval = 1
           are set by the OpenFOAM coupling: the labels are transported every solver time step */
    }

    void user_set_phase_dict(void *handle)
    {
        rcfd_set_dict_property_i(handle, "Phase_Dict", -1, ph_gas, -1, "number_of_data", 1);

        rcfd_set_dict_property_r(handle, "Phase_Dict", -1, ph_gas, -1, "density", 1.2);

        rcfd_set_dict_property_i(handle, "Phase_Dict", -1, ph_solid, -1, "number_of_data", 2);

        rcfd_set_dict_property_r(handle, "Phase_Dict", -1, ph_solid, -1, "density", 2500.0);
    }

    void user_set_rec_dict(void *handle)
    {
        rcfd_set_dict_property_s(handle, "Rec_Dict", -1, -1, -1, "format", "replay");
    }

    void user_set_cell_dict(void *handle, const short i_layer)
    {
        fluent_set_cell_dict(handle, i_layer);      /* here: cell counts of the OpenFOAM mesh */

        rcfd_set_dict_property_r(handle, "Cell_Dict", i_layer, -1, -1, "next_neigh_cell_skin_factor", 1.5);
    }

    void user_set_face_dict(void *handle, const short i_layer)
    {
        fluent_set_face_dict(handle, i_layer);      /* here: face counts of the OpenFOAM mesh */
    }

    /* register user functions */
    /*************************************************************************************/

    void user_register_functions(void *handle)
    {
        REGISTER_RCFD_UDF( user_set_solver_dict );
        REGISTER_RCFD_UDF( user_set_phase_dict );
        REGISTER_RCFD_UDF( user_set_rec_dict );

        REGISTER_RCFD_UDF( user_set_cell_dict );
        REGISTER_RCFD_UDF( user_set_face_dict );
    }

#endif
