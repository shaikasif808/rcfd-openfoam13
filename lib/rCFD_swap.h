/******************************************************************************
 * @internal
 * @file        rCFD_swap.h
 * @copyright   (C) 2023
 *              Particulate Flow Modelling,
 *              Johannes Kepler University, Linz, Austria
 * @website     www.particulate-flow.at
 * @authors     Stefan Pirker
 *****************************************************************************/

#ifndef RCFD_SWAP
#define RCFD_SWAP

#include "rCFD_types.h"
#include "rCFD_macros.h"
#include "rCFD_topo.h"
#include "rCFD_parallel.h"

/* TUHH 11/2026: refresh the corona (ghost) copies of the phase data between drift
 * sub-loops (number_of_drift_loops > 1). 0 = off (ghost copies are updated by the
 * communicated exchange only), 1 = on (parallel_map_phase_data before every sub-loop
 * after the first). Diagnostic: a k = 0 run must give the same result for 1 and N loops. */
#ifndef RCFD_DRIFT_LOOP_REFRESH
#define RCFD_DRIFT_LOOP_REFRESH 0
#endif

/* TUHH 11/2026: per-sub-loop diagnostics of the conservative drift, printed on node 0
 * every RCFD_DRIFT_DIAG_EVERY time steps (0 = off): summed |exchange| per sub-loop,
 * number of faces touched by the limiter, net flux of data 0 (positive = c0 -> c1). */
#ifndef RCFD_DRIFT_DIAG_EVERY
#define RCFD_DRIFT_DIAG_EVERY 40
#endif

    void swap_init(rCFD_type* rcfd, const int i_phase)
    {
#if RP_NODE
        int i_layer = rcfd->Solver->current_layer;

        int cell_count = _CELL_DICT(i_layer).number_of_cells;
        int data_count = rcfd->Phase_Dict[i_phase].number_of_data;

        memset(_C(i_layer).weight_after_swap, 0, cell_count *              sizeof(double));
        memset(_C(i_layer).data_swap[0],      0, cell_count * data_count * sizeof(double));

        parallel_map_phase_data(rcfd, i_phase);
#endif
    }

    void swap_data(rCFD_type* rcfd, const int i_phase)
    {
#if RP_NODE

#define LOCAL_DEBUG_FLAG    0

#if 1   // local vars
    int     i_layer, i_data, i_while, i_diff, i_frame, i_frame_c0, i_frame_c1, i_face, i_cell;
    int     number_of_diff_loops, c, c0, c1, swap_exchange_active;
    int     *active_swap_face = NULL;
    double  swap_per_loop, data_exchange, average_local_diff, local_cell_distance, local_face_area, max_local_swap_exchange;
    double  vol_c0, vol_c1, vol_flip, max_vol_flip;
    double  *local_swap_exchange = NULL;
#endif

    i_layer = rcfd->Solver->current_layer;
    i_frame = rcfd->Solver->current_frame;

    i_frame_c0 = i_frame;
    i_frame_c1 = i_frame;

    swap_init(rcfd, i_phase);

    if(rcfd->Solver_Dict->local_diff_on){

        for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

            switch (rcfd->Data_Dict[i_phase][i_data].type){

                case RCFD_DATA_CONCENTRATION:
                case RCFD_DATA_TEMPERATURE:
                {
                    // init swap_exchange_active, calc swap exchange per face[i_face], set active_swap_face[i_face]
                    {
                        // while loop
                        swap_exchange_active = 1;

                        i_while = 0;

                        // swap_exchange_per_face, active_swap_face
                        local_swap_exchange = (double*)malloc(_FACE_DICT(i_layer).number_of_faces * sizeof(double));

                        active_swap_face = (int*)malloc(_FACE_DICT(i_layer).number_of_faces * sizeof(int));

                        max_local_swap_exchange = 0.0;

                        for(i_face = 0; i_face < _FACE_DICT(i_layer).number_of_faces; i_face++){

                            c0 = _F(i_layer).c0[i_face];
                            c1 = _F(i_layer).c1[i_face];

                            /* TUHH 11/2026: boundary face to a corona cell without owner, see rCFD_parallel.h */
                            if(parallel_ext_face_blocked(i_face)){

                                local_swap_exchange[i_face] = 0.0;
                                active_swap_face[i_face] = 0;

                                continue;
                            }

                            average_local_diff = (_C(i_layer).local_diff[i_frame_c0][i_phase][c0] + _C(i_layer).local_diff[i_frame_c1][i_phase][c1]) / 2.0;

                            local_cell_distance = cell2cell_distance(rcfd, i_layer, c0, c1);

                            local_face_area = topo_face_area(rcfd, i_layer, i_face);

                            local_swap_exchange[i_face] =     average_local_diff
                                                            * local_face_area
                                                            / local_cell_distance
                                                            * rcfd->Solver_Dict->global_time_step;   // (m3)

                            if(local_swap_exchange[i_face] > max_local_swap_exchange){

                                max_local_swap_exchange = local_swap_exchange[i_face];
                            }

                            if(fabs(local_swap_exchange[i_face]) > 0.0){

                                active_swap_face[i_face] = 1;
                            }
                            else{

                                active_swap_face[i_face] = 0;
                            }
                        }

                        max_local_swap_exchange = PRF_GRSUM1(max_local_swap_exchange);

                        if(LOCAL_DEBUG_FLAG){

                            report_debug0(rcfd, RCFD_REPORT_CATEGORY_RUN_SWAP,
                                          "swap_data: max_local_swap_exchange = %e\n",
                                          max_local_swap_exchange);
                        }
                    }

                    while((swap_exchange_active) && (i_while < rcfd->Solver_Dict->face_swap_max_loops)){

                        // init loop vars, _C(i_layer).data_swap
                        {
                            i_while++;

                            swap_exchange_active = 0;

                            for(i_cell = 0; i_cell < _CELL_DICT(i_layer).number_of_cells; i_cell++){

                                _C(i_layer).data_swap[i_cell][i_data] = 0.0;
                            }
                        }

                        for(i_face = 0; i_face < _FACE_DICT(i_layer).number_of_faces; i_face++){

                            if(active_swap_face[i_face]){

                                // get neighboring cells, vof's, max vol_flip
                                {
                                    c0 = _F(i_layer).c0[i_face];
                                    c1 = _F(i_layer).c1[i_face];

                                    vol_c0 = _C(i_layer).volume[c0] * _C(i_layer).vof[i_frame_c0][c0][i_phase];
                                    vol_c1 = _C(i_layer).volume[c1] * _C(i_layer).vof[i_frame_c1][c1][i_phase];

                                    if(vol_c0 < vol_c1){

                                        max_vol_flip = vol_c0 * rcfd->Solver_Dict->face_swap_max_per_loop;
                                    }
                                    else{

                                        max_vol_flip = vol_c1 * rcfd->Solver_Dict->face_swap_max_per_loop;
                                    }
                                }

                                // calc. vol_flip
                                {
                                    if(local_swap_exchange[i_face] <= max_vol_flip){

                                        vol_flip = local_swap_exchange[i_face];
                                        local_swap_exchange[i_face] = 0.0;
                                        active_swap_face[i_face] = 0;
                                    }
                                    else{

                                        vol_flip = max_vol_flip;
                                        local_swap_exchange[i_face] -= max_vol_flip;
                                    }
                                }

                                // calc. data_exchange
                                data_exchange = (_C(i_layer).data[i_phase][c1][i_data] - _C(i_layer).data[i_phase][c0][i_data]) / 2.;

                                // update _C(i_layer).data_swap
                                {
                                    _C(i_layer).data_swap[c1][i_data] -= data_exchange * vol_flip;
                                    _C(i_layer).data_swap[c0][i_data] += data_exchange * vol_flip;
                                }

                                // update active_swap_face, swap_exchange_active
                                {
                                    if(local_swap_exchange[i_face] > 0.0){

                                        swap_exchange_active++;
                                    }
                                    else{

                                        active_swap_face[i_face] = 0;
                                    }
                                }
                            }   // if active face
                        }

                        // update _C(i_layer).data
                        {
                            for(i_cell = 0; i_cell < _CELL_DICT(i_layer).number_of_cells; i_cell++){

                            if((_C(i_layer).volume[i_cell] * _C(i_layer).vof[i_frame][i_cell][i_phase]) > 0.0){

                                    _C(i_layer).data[i_phase][i_cell][i_data] += _C(i_layer).data_swap[i_cell][i_data] / (_C(i_layer).volume[i_cell] * _C(i_layer).vof[i_frame][i_cell][i_phase]);
                                }
                            }
                        }

                        // map host data to corona cells
                        if(rcfd->Solver_Dict->parallel_diffusion_on && i_layer == 0){

                            parallel_communicate_data(rcfd, i_layer, i_phase, i_data, RCFD_PARALLEL_MAP_HOST_FORMAT);
                        }

                        swap_exchange_active = PRF_GISUM1(swap_exchange_active);

                        // local debug
                        if(LOCAL_DEBUG_FLAG){

                            report_debug0(rcfd, RCFD_REPORT_CATEGORY_RUN_SWAP,
                                          "data_swap: i_while = %d, swap_exchange_active = %d",
                                          i_while, swap_exchange_active);
                        }
                    }   // while loop

                    // free local vars
                    {
                        free( local_swap_exchange);
                        free( active_swap_face);
                    }

                    break;
                }
            }
        }
    }
    else{   // global diff

        for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

            switch (rcfd->Data_Dict[i_phase][i_data].type){

                case RCFD_DATA_CONCENTRATION:
                case RCFD_DATA_TEMPERATURE:
                {

                    // define global number_of_diff_loops, swap_per_loop
                    {
                        if(fabs(rcfd->Data_Dict[i_phase][i_data].face_swap_diffusion) > rcfd->Solver_Dict->face_swap_max_per_loop){

                            number_of_diff_loops = (int)(fabs(rcfd->Data_Dict[i_phase][i_data].face_swap_diffusion) / rcfd->Solver_Dict->face_swap_max_per_loop);  // ( > 1)
                            swap_per_loop = rcfd->Data_Dict[i_phase][i_data].face_swap_diffusion / (double)number_of_diff_loops;
                        }
                        else{

                            number_of_diff_loops = 1;
                            swap_per_loop = rcfd->Data_Dict[i_phase][i_data].face_swap_diffusion;
                        }
                    }

                    // perform face swaps until global diffusive fluxes are achieved
                    for(i_diff = 0; i_diff < number_of_diff_loops; i_diff++){

                        /* init. data_swap */
                        for(i_cell = 0; i_cell < _CELL_DICT(i_layer).number_of_int_cells; i_cell++){

                            _C(i_layer).data_swap[i_cell][i_data] = 0.0;
                        }

                        /* define swap masses */
                        for(i_face = 0; i_face < _FACE_DICT(i_layer).number_of_faces; i_face++){

                            /* TUHH 11/2026: boundary face to a corona cell without owner, see rCFD_parallel.h */
                            if(parallel_ext_face_blocked(i_face)) continue;

                            c0 = _F(i_layer).c0[i_face];
                            c1 = _F(i_layer).c1[i_face];

                            if(_C(i_layer).data[i_phase][c0][i_data] > _C(i_layer).data[i_phase][c1][i_data]){

                                c = c0; c0 = c1; c1 = c;
                            }

                            if(_C(i_layer).data[i_phase][c1][i_data] > _C(i_layer).data[i_phase][c0][i_data]){

                                // flip volume
                                {
                                    if((_C(i_layer).volume[c0] * _C(i_layer).vof[i_frame_c0][c0][i_phase]) < (_C(i_layer).volume[c1] * _C(i_layer).vof[i_frame_c1][c1][i_phase])){

                                        vol_flip = _C(i_layer).volume[c0] * _C(i_layer).vof[i_frame_c0][c0][i_phase];
                                    }
                                    else{

                                        vol_flip = _C(i_layer).volume[c1] * _C(i_layer).vof[i_frame_c1][c1][i_phase];
                                    }
                                }

                                // data_exchange
                                {
                                    data_exchange = (_C(i_layer).data[i_phase][c1][i_data] - _C(i_layer).data[i_phase][c0][i_data]) / 2.;

                                    if(swap_per_loop < 0.0){    // neg. diff.

                                        if((1.0 - _C(i_layer).data[i_phase][c1][i_data]) < _C(i_layer).data[i_phase][c0][i_data]){

                                            data_exchange = 1.0 - _C(i_layer).data[i_phase][c1][i_data];
                                        }
                                        else{

                                            data_exchange = _C(i_layer).data[i_phase][c0][i_data];
                                        }
                                    }
                                }

                                // _C(i_layer).data_swap
                                {
                                    _C(i_layer).data_swap[c1][i_data] -= vol_flip * data_exchange * swap_per_loop;
                                    _C(i_layer).data_swap[c0][i_data] += vol_flip * data_exchange * swap_per_loop;
                                }
                            }
                        }

                        /* update data by data_swap */
                        for(i_cell = 0; i_cell < _CELL_DICT(i_layer).number_of_int_cells; i_cell++){

                            if((_C(i_layer).volume[i_cell] * _C(i_layer).vof[i_frame][i_cell][i_phase]) > 0.0){

                                _C(i_layer).data[i_phase][i_cell][i_data] += _C(i_layer).data_swap[i_cell][i_data] / (_C(i_layer).volume[i_cell] * _C(i_layer).vof[i_frame][i_cell][i_phase]);
                            }
                        }

                        /* parallel diffusion: map host data to other nodes' corona cells */
                        if(rcfd->Solver_Dict->parallel_diffusion_on && i_layer == 0){

                            parallel_communicate_data(rcfd, i_layer, i_phase, i_data, RCFD_PARALLEL_MAP_HOST_FORMAT);
                        }
                    }

                    break;
                }
            }
        }
    }

#undef  LOCAL_DEBUG_FLAG

#endif
    }

    void swap_drift_data(rCFD_type* rcfd, const int i_phase)
    {
#if RP_NODE

    /* Existing Limitations (11/21)
        L1: drifting pattern might depend on face order
        L2: drifting of tmperature_data not yet implemented
        L3: drifting only implemented for i_layer == 0
        L4: local drift only implemented for concentration_data
    */

#define LOCAL_DEBUG_FLAG    0

#if 1   // local vars
        int     i_layer, i_data, i_face, i_drift, i_cell, i_dim, i_frame;
        int     c0, c1, c, i_frame_c0, i_frame_c1, current_data_index;
        double  d32_fraction, d32_bulk, d32_local_nom, d32_local_denom, d32_local, local_conc;
        double  local_mass_c0, local_mass_c1, local_mass, drift_volume, hindering_factor;
        double  current_conc, sum_of_other_concs, conc_error;
        double  local_drift_exchange, available_c1_mass;
        int     int_cell_count, drift_face_count;   /* TUHH 11/2026, MPI fix */
        double  vof_cut;
#endif

        i_layer = rcfd->Solver->current_layer;
        i_frame = rcfd->Solver->current_frame;

        i_frame_c0 = i_frame;
        i_frame_c1 = i_frame;

        swap_init(rcfd, i_phase);

        /* ------------------------------------------------------------------
         * TUHH 11/2026: drift across partition boundaries (MPI fix).
         *
         * A partition-boundary face is in the face list of BOTH neighbouring
         * nodes (c0 internal, c1 ext/corona cell). The former code accumulated
         * the face exchange into both cells on both nodes and then called
         * parallel_communicate(SUM_CELLS), which adds the corona values to the
         * hosting cells: every boundary flux reached each owner cell twice.
         * It also relied on corona-cell vof (never filled in run mode, see
         * parallel_expand_initial_fields) and on the orientation of the ext
         * face areas (not flipped with c0/c1, see topo_orient_face_areas).
         *
         * Now, like swap_data does for the face swaps: every node evaluates
         * the flux of each of its faces from identical inputs (data, vof and
         * volume of corona cells are copies of the owner values), applies it
         * to its INTERNAL cells only, and the corona copies of the data are
         * refreshed afterwards. Each flux is applied exactly once per owner
         * cell, species mass is conserved across partitions and the drift
         * step does not depend on the partitioning. No SUM_CELLS communication.
         * (The balance correction in rCFD_balance.h still works on internal
         * faces only; whatever mass error it has to repair is repaired in a
         * partition-dependent way. That is not changed here.)
         *
         * parallel_diffusion_on == 0 now means: no drift across partition
         * boundaries (ext faces are skipped).
         *
         * A face exchanges only if BOTH cells are above vof_cutoff_limit, the
         * same threshold the cell update uses; before, a cell at 0 < vof <=
         * cutoff took part in the face exchange but was never updated, which
         * created/destroyed species mass at the bed surface.
         * ------------------------------------------------------------------ */
        int_cell_count   = _CELL_DICT(i_layer).number_of_int_cells;
        drift_face_count = (rcfd->Solver_Dict->parallel_diffusion_on) ? _FACE_DICT(i_layer).number_of_faces
                                                                      : _FACE_DICT(i_layer).number_of_int_faces;
        vof_cut          = rcfd->Solver_Dict->vof_cutoff_limit;

        // set local_drift_adaption per face
        {
            /* Phase_Dict.local_drift_on selects how the per-face multiplier of
             * Data_Dict.drift_velocity is built from the LOCAL Sauter diameter
             * d32_local of the two cells sharing the face:
             *
             *   1 : original rCFD-PFM behaviour. drift_velocity[i] is the drift of
             *       fraction i evaluated at the bulk reference Phase_Dict.d32_diam;
             *       the multiplier rescales it linearly to the local reference,
             *           (d_i - d32_local) / (d_i - d32_bulk).
             *
             *   2 : local square law (TUHH C04, 11/2026). drift_velocity holds
             *       k_d * e, the SAME vector for every fraction (k_d in 1/(m s),
             *       e a direction vector). The multiplier is
             *           (d_i - d32_local) * |d_i - d32_local|      [m^2]
             *       so the applied drift is u_i = k_d e (d_i - d32_local)|d_i - d32_local|.
             *       With e pointing against gravity-sense, i.e. DOWN, fractions coarser
             *       than the local Sauter diameter sink and finer ones rise, and a
             *       locally monodisperse mixture does not drift at all. d32_diam is
             *       not used. This is the form needed when the PSD evolves during the
             *       run (granulation), because no fixed reference exists.
             *
             *   3 : local linear law, multiplier (d_i - d32_local)            [m].
             *
             * d32_local: concentration data are MASS fractions, hence
             *       1 / d32 = sum_i ( c_i / d_i ).
             * The earlier sum(c_i d_i^3) / sum(c_i d_i^2) is the NUMBER-fraction form
             * and over-weights the coarse fractions when fed with mass fractions.
             * A face without any concentration information gets multiplier 0 (no
             * drift) instead of reading Data_Dict out of range. */

            if(rcfd->Phase_Dict[i_phase].local_drift_on){

                d32_bulk = rcfd->Phase_Dict[i_phase].d32_diam;

                for(i_face = 0; i_face < _FACE_DICT(i_layer).number_of_faces; i_face++){

                    c0 = _F(i_layer).c0[i_face];
                    c1 = _F(i_layer).c1[i_face];

                    d32_local_nom   = 0.0;      /* sum c_i        */
                    d32_local_denom = 0.0;      /* sum c_i / d_i  */

                    for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

                        if((rcfd->Data_Dict[i_phase][i_data].type == RCFD_DATA_CONCENTRATION) &&
                           (rcfd->Data_Dict[i_phase][i_data].d1_diam > 0.0)){

                            local_conc = (_C(i_layer).data[i_phase][c0][i_data] + _C(i_layer).data[i_phase][c1][i_data]) / 2.0;

                            if(local_conc > 0.0){

                                d32_local_nom   += local_conc;
                                d32_local_denom += local_conc / rcfd->Data_Dict[i_phase][i_data].d1_diam;
                            }
                        }
                    }

                    if(d32_local_denom > 0.0){

                        d32_local = d32_local_nom / d32_local_denom;
                    }
                    else{

                        d32_local = -1.0;   /* flag: no PSD information on this face */
                    }

                    for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

                        d32_fraction = rcfd->Data_Dict[i_phase][i_data].d1_diam;

                        if(d32_local <= 0.0){

                            _F(i_layer).local_drift_adaption[i_face][i_data] = 0.0;
                        }
                        else if(rcfd->Phase_Dict[i_phase].local_drift_on == 2){

                            _F(i_layer).local_drift_adaption[i_face][i_data] =   (d32_fraction - d32_local)
                                                                               * fabs(d32_fraction - d32_local);
                        }
                        else if(rcfd->Phase_Dict[i_phase].local_drift_on == 3){

                            _F(i_layer).local_drift_adaption[i_face][i_data] = (d32_fraction - d32_local);
                        }
                        else{   /* == 1, original form */

                            if(d32_fraction != d32_bulk){

                                _F(i_layer).local_drift_adaption[i_face][i_data] = (d32_fraction - d32_local) / (d32_fraction - d32_bulk);
                            }
                            else{

                                _F(i_layer).local_drift_adaption[i_face][i_data] = 1.0;
                            }
                        }
                    }
                }
            }
        }

        /* ------------------------------------------------------------------
         * TUHH 11/2026: conservative multi-species drift for concentration data
         * (Phase_Dict.conservative_drift_on == 1: donor vof, == 2: min vof at the face).
         *
         * All drifting concentration data of the phase are exchanged together.
         * Per face, the species volume fluxes phi_i = c_i(donor) vof(donor) u_i.A dt
         * are reduced by their mixture value, F_i = phi_i - c_i(mixture donor) sum_j phi_j,
         * so that sum_i F_i == 0 exactly: what one species gains at a face the
         * others lose, species mass is conserved, the concentrations keep summing
         * to one without any re-normalisation, and no centre-of-mass correction
         * is needed (only velocity differences act). The original sequential
         * per-data exchange with drift_control_conc_sum_on creates/destroys
         * species mass of the order u dt / dx per step, which the balance
         * correction then redistributes along the concentration gradients.
         *
         * hindered_drift_on == 2 scales the face drift by
         * min(1, vof_face / drift_dilute_vof): no size drift in the dilute
         * splash zone, where particles fly ballistically.
         * ------------------------------------------------------------------ */
        if(rcfd->Phase_Dict[i_phase].conservative_drift_on){

            int     n_data     = rcfd->Phase_Dict[i_phase].number_of_data;
            int     cell_count = _CELL_DICT(i_layer).number_of_cells;
            int     n_loops    = rcfd->Solver_Dict->number_of_drift_loops;
            int     any_drift  = 0, donor_mix, k_data;
            double  *vol_i     = (double*)malloc(n_data * sizeof(double));
            double  *flux_i    = (double*)malloc(n_data * sizeof(double));
            double  rho        = rcfd->Phase_Dict[i_phase].density;
            double  vof0, vof1, vol0, vol1, vof_face, vof_x, dilute, vol_mix, limit, scale, local_m, sum_c;
            double  diag_abs = 0.0, diag_net0 = 0.0, diag_lim = 0.0, diag_faces = 0.0;   /* per sub-loop, this node */
            double  diag_bnd_net = 0.0, diag_bnd_abs = 0.0;

            for(i_data = 0; i_data < n_data; i_data++){

                if(rcfd->Data_Dict[i_phase][i_data].drifting_data &&
                   (rcfd->Data_Dict[i_phase][i_data].type == RCFD_DATA_CONCENTRATION)){

                    any_drift = 1;
                }
            }

            if(any_drift && vol_i && flux_i){

                for(i_drift = 0; i_drift < n_loops; i_drift++){

                    if(RCFD_DRIFT_LOOP_REFRESH && (i_drift > 0)){

                        parallel_map_phase_data(rcfd, i_phase);
                    }

                    /* data_swap (free at this point of the time step) accumulates the
                     * species mass exchange [kg] per cell */
                    for(i_cell = 0; i_cell < cell_count; i_cell++){

                        for(i_data = 0; i_data < n_data; i_data++){

                            _C(i_layer).data_swap[i_cell][i_data] = 0.0;
                        }
                    }

                    diag_abs = 0.0; diag_net0 = 0.0; diag_lim = 0.0; diag_faces = 0.0;
                    diag_bnd_net = 0.0; diag_bnd_abs = 0.0;

                    for(i_face = 0; i_face < drift_face_count; i_face++){

                        double face_w;  /* diagnostics only: an ext face is seen by two nodes */

                        if(parallel_ext_face_blocked(i_face)) continue;

                        c0 = _F(i_layer).c0[i_face];
                        c1 = _F(i_layer).c1[i_face];

                        face_w = (c1 < int_cell_count) ? 1.0 : 0.5;

                        vof0 = _C(i_layer).vof[i_frame][c0][i_phase];
                        vof1 = _C(i_layer).vof[i_frame][c1][i_phase];
                        vol0 = _C(i_layer).volume[c0];
                        vol1 = _C(i_layer).volume[c1];

                        /* don't drift across interfaces; both cells must be cells the update acts on */
                        if((vol0 <= 0.0) || (vol1 <= 0.0) || (vof0 <= vof_cut) || (vof1 <= vof_cut)) continue;

                        vof_face = 0.5 * (vof0 + vof1);

                        dilute = 1.0;

                        if(rcfd->Phase_Dict[i_phase].hindered_drift_on == 1){

                            dilute = 1.0 - vof_face / rcfd->Phase_Dict[i_phase].vof_max;

                            if(dilute < 0.0) dilute = 0.0;
                        }
                        else if((rcfd->Phase_Dict[i_phase].hindered_drift_on == 2) &&
                                (rcfd->Phase_Dict[i_phase].drift_dilute_vof > 0.0)){

                            dilute = vof_face / rcfd->Phase_Dict[i_phase].drift_dilute_vof;

                            if(dilute > 1.0) dilute = 1.0;
                        }

                        if(dilute <= 0.0) continue;

                        /* species volume fluxes c0 -> c1 (signed), donor-cell concentration and vof */
                        vol_mix = 0.0;

                        for(i_data = 0; i_data < n_data; i_data++){

                            vol_i[i_data]  = 0.0;
                            flux_i[i_data] = 0.0;

                            if(!(rcfd->Data_Dict[i_phase][i_data].drifting_data &&
                                 (rcfd->Data_Dict[i_phase][i_data].type == RCFD_DATA_CONCENTRATION))) continue;

                            for(i_dim = 0; i_dim < 3; i_dim++){

                                vol_i[i_data] +=   _F(i_layer).local_drift_adaption[i_face][i_data]
                                                 * rcfd->Data_Dict[i_phase][i_data].drift_velocity[i_dim]
                                                 * _F(i_layer).area[i_face][i_dim]
                                                 * rcfd->Solver->timestep_width[i_layer]
                                                 / (double)n_loops;
                            }

                            vol_i[i_data] *= dilute;

                            /* conservative_drift_on == 1: donor-cell vof scales the flux;
                             * conservative_drift_on == 2: the smaller vof of the two cells does, so a
                             *   drift exchange across the bed surface (dense below, dilute above) is
                             *   limited by what the dilute side can supply instead of draining it. */
                            if(rcfd->Phase_Dict[i_phase].conservative_drift_on == 2){

                                vof_x = (vof0 < vof1) ? vof0 : vof1;
                            }
                            else{

                                vof_x = (vol_i[i_data] > 0.0) ? vof0 : vof1;
                            }

                            if(vol_i[i_data] > 0.0){

                                flux_i[i_data] = _C(i_layer).data[i_phase][c0][i_data] * vof_x * vol_i[i_data];
                            }
                            else{

                                flux_i[i_data] = _C(i_layer).data[i_phase][c1][i_data] * vof_x * vol_i[i_data];
                            }

                            vol_mix += flux_i[i_data];
                        }

                        /* mixture-relative fluxes: sum_i F_i = 0 because the donor concentrations sum to 1 */
                        donor_mix = (vol_mix > 0.0) ? c0 : c1;

                        for(i_data = 0; i_data < n_data; i_data++){

                            if(!(rcfd->Data_Dict[i_phase][i_data].drifting_data &&
                                 (rcfd->Data_Dict[i_phase][i_data].type == RCFD_DATA_CONCENTRATION))) continue;

                            flux_i[i_data] -= _C(i_layer).data[i_phase][donor_mix][i_data] * vol_mix;
                        }

                        /* limiter: through one face a species loses at most half of its content
                         * in the giving cell; the whole face vector is scaled to keep sum_i F_i = 0 */
                        scale = 1.0;

                        for(i_data = 0; i_data < n_data; i_data++){

                            if(flux_i[i_data] > 0.0){

                                limit = 0.5 * _C(i_layer).data[i_phase][c0][i_data] * vof0 * vol0;

                                if(flux_i[i_data] > limit){

                                    if(limit / flux_i[i_data] < scale) scale = limit / flux_i[i_data];
                                }
                            }
                            else if(flux_i[i_data] < 0.0){

                                limit = 0.5 * _C(i_layer).data[i_phase][c1][i_data] * vof1 * vol1;

                                if(-flux_i[i_data] > limit){

                                    if(limit / (-flux_i[i_data]) < scale) scale = limit / (-flux_i[i_data]);
                                }
                            }
                        }

                        if(scale < 0.0) scale = 0.0;

                        if(scale < 1.0) diag_lim += face_w;
                        diag_faces += face_w;

                        for(i_data = 0; i_data < n_data; i_data++){

                            if(flux_i[i_data] != 0.0){

                                /* internal (owner) cells only; the neighbour node applies the other side */
                                if(c0 < int_cell_count) _C(i_layer).data_swap[c0][i_data] -= flux_i[i_data] * scale * rho;
                                if(c1 < int_cell_count) _C(i_layer).data_swap[c1][i_data] += flux_i[i_data] * scale * rho;

                                diag_abs += face_w * fabs(flux_i[i_data] * scale * rho);
                                if(i_data == 0) diag_net0 += face_w * flux_i[i_data] * scale * rho;

                                /* partition-boundary balance: what leaves the internal cell here must
                                 * enter the internal cell of the neighbour node; the global sum is zero */
                                if(c1 >= int_cell_count){

                                    diag_bnd_net += flux_i[i_data] * scale * rho;
                                    diag_bnd_abs += 0.5 * fabs(flux_i[i_data] * scale * rho);
                                }
                            }
                        }
                    }

                    if(RCFD_DRIFT_DIAG_EVERY > 0 && (rcfd->Solver->global_run_counter % RCFD_DRIFT_DIAG_EVERY) == 0){

                        double g_abs = PRF_GRSUM1(diag_abs), g_net0 = PRF_GRSUM1(diag_net0);
                        double g_lim = PRF_GRSUM1(diag_lim), g_faces = PRF_GRSUM1(diag_faces);

                        double g_bnd_net = PRF_GRSUM1(diag_bnd_net), g_bnd_abs = PRF_GRSUM1(diag_bnd_abs);

                        report_print0(rcfd, "[drift] step %d phase %d loop %d/%d: sum|dm| = %.4e kg, net data0 (c0->c1) = %+.4e kg, limited faces %.0f of %.0f",
                                      rcfd->Solver->global_run_counter, i_phase, i_drift + 1, n_loops, g_abs, g_net0, g_lim, g_faces);

                        report_print0(rcfd, "[drift] step %d phase %d loop %d/%d: partition boundaries exchange %.4e kg, imbalance %+.3e kg (must be round-off)",
                                      rcfd->Solver->global_run_counter, i_phase, i_drift + 1, n_loops, g_bnd_abs, g_bnd_net);
                    }

                    /* no SUM_CELLS exchange (it double counted the partition-boundary fluxes):
                     * data_swap of the internal cells is complete, see the note at the top */

                    /* apply to internal cells; sum over species of the exchange is zero per
                     * cell, so the concentrations stay normalised up to round-off */
                    for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                        if(_C(i_layer).vof[i_frame][i_cell][i_phase] > rcfd->Solver_Dict->vof_cutoff_limit){

                            local_m = _C(i_layer).volume[i_cell] * _C(i_layer).vof[i_frame][i_cell][i_phase] * rho;

                            sum_c = 0.0;

                            for(i_data = 0; i_data < n_data; i_data++){

                                if(rcfd->Data_Dict[i_phase][i_data].type == RCFD_DATA_CONCENTRATION){

                                    _C(i_layer).data[i_phase][i_cell][i_data] += _C(i_layer).data_swap[i_cell][i_data] / local_m;

                                    if(_C(i_layer).data[i_phase][i_cell][i_data] < 0.0) _C(i_layer).data[i_phase][i_cell][i_data] = 0.0;

                                    sum_c += _C(i_layer).data[i_phase][i_cell][i_data];
                                }
                            }

                            if(sum_c > 1.0e-12){

                                for(k_data = 0; k_data < n_data; k_data++){

                                    if(rcfd->Data_Dict[i_phase][k_data].type == RCFD_DATA_CONCENTRATION){

                                        _C(i_layer).data[i_phase][i_cell][k_data] /= sum_c;
                                    }
                                }
                            }
                        }
                    }

                    /* corona copies follow their owner cells (next sub-loop, sequential part) */
                    if(rcfd->Solver_Dict->parallel_diffusion_on){

                        parallel_map_phase_data(rcfd, i_phase);
                    }
                }
            }

            if(vol_i)  free(vol_i);
            if(flux_i) free(flux_i);
        }

        // for each data calculate drift exchange
        for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

            /* concentration data are already done above in conservative mode */
            if(rcfd->Phase_Dict[i_phase].conservative_drift_on &&
               (rcfd->Data_Dict[i_phase][i_data].type == RCFD_DATA_CONCENTRATION)) continue;

            if(rcfd->Data_Dict[i_phase][i_data].drifting_data){

                for(i_drift = 0; i_drift < rcfd->Solver_Dict->number_of_drift_loops; i_drift++){

                    /* TUHH 11/2026: refresh the corona (external) cell copies of
                     * the phase data before every sub-loop after the first. The
                     * first sub-loop works on the copies made by swap_init; later
                     * sub-loops would otherwise read stale values at partition
                     * boundaries. With number_of_drift_loops == 1 nothing changes. */
                    if(RCFD_DRIFT_LOOP_REFRESH && (i_drift > 0)){

                        parallel_map_phase_data(rcfd, i_phase);
                    }

                    /* init _C(i_layer).drift_exchange */
                    for(i_cell = 0; i_cell < _CELL_DICT(i_layer).number_of_cells; i_cell++){

                        _C(i_layer).drift_exchange[i_cell] = 0.0;
                    }

                    switch (rcfd->Data_Dict[i_phase][i_data].type){

                        case RCFD_DATA_GENERIC:
                        {
                            /* set _C(i_layer).drift_exchange */
                            for(i_face = 0; i_face < drift_face_count; i_face++){

                                if(parallel_ext_face_blocked(i_face)) continue;

                                c0 = _F(i_layer).c0[i_face];
                                c1 = _F(i_layer).c1[i_face];

                                local_mass_c0 = _C(i_layer).volume[c0] * _C(i_layer).vof[i_frame_c0][c0][i_phase] * rcfd->Phase_Dict[i_phase].density;
                                local_mass_c1 = _C(i_layer).volume[c1] * _C(i_layer).vof[i_frame_c1][c1][i_phase] * rcfd->Phase_Dict[i_phase].density;

                                /* don't drift across interfaces */
                                if(local_mass_c0 * local_mass_c1 > 0.0){

                                    drift_volume = 0.0;

                                    for(i_dim = 0; i_dim < 3; i_dim++){

                                        drift_volume +=   rcfd->Data_Dict[i_phase][i_data].drift_velocity[i_dim]
                                                        * _F(i_layer).area[i_face][i_dim]
                                                        * rcfd->Solver->timestep_width[i_layer]
                                                        / (double)rcfd->Solver_Dict->number_of_drift_loops;
                                    }

                                    /* no flux through this face (also keeps both nodes of an ext face on the same donor) */
                                    if(drift_volume == 0.0) continue;

                                    /* flip cells, such that flux is from c0 to c1 */
                                    if(drift_volume < 0.0){

                                        c = c0; c0 = c1; c1 = c;

                                        i_frame = i_frame_c0; i_frame_c0 = i_frame_c1; i_frame_c1 = i_frame;

                                        local_mass = local_mass_c0; local_mass_c0 = local_mass_c1; local_mass_c1 = local_mass;

                                        drift_volume *= -1.0;
                                    }

                                    if(drift_volume > _C(i_layer).volume[c0]){

                                        drift_volume = _C(i_layer).volume[c0];
                                    }

                                    local_drift_exchange = 0.0;

                                    if(_C(i_layer).volume[c0] > 0.0){

                                        local_drift_exchange = _C(i_layer).data[i_phase][c0][i_data] * drift_volume /_C(i_layer).volume[c0];
                                    }

                                    /* internal (owner) cells only, see the note at the top */
                                    if(c0 < int_cell_count) _C(i_layer).drift_exchange[c0] -= local_drift_exchange;
                                    if(c1 < int_cell_count) _C(i_layer).drift_exchange[c1] += local_drift_exchange;
                                }
                            }

                            /* adjust data of internal cells by _C(i_layer).drift_exchange */
                            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                                _C(i_layer).data[i_phase][i_cell][i_data] += _C(i_layer).drift_exchange[i_cell];
                            }

                            /* corona copies follow their owner cells */
                            if(rcfd->Solver_Dict->parallel_diffusion_on){

                                parallel_communicate_data(rcfd, i_layer, i_phase, i_data, RCFD_PARALLEL_MAP_HOST_FORMAT);
                            }

                            break;
                        }
                        case RCFD_DATA_CONCENTRATION:
                        {
                            double seq_bnd_net = 0.0, seq_bnd_abs = 0.0;   /* diagnostics, this node */

                            /* set _C(i_layer).drift_exchange */
                            for(i_face = 0; i_face < drift_face_count; i_face++){

                                if(parallel_ext_face_blocked(i_face)) continue;

                                c0 = _F(i_layer).c0[i_face];
                                c1 = _F(i_layer).c1[i_face];

                                local_mass_c0 = _C(i_layer).volume[c0] * _C(i_layer).vof[i_frame_c0][c0][i_phase] * rcfd->Phase_Dict[i_phase].density;
                                local_mass_c1 = _C(i_layer).volume[c1] * _C(i_layer).vof[i_frame_c1][c1][i_phase] * rcfd->Phase_Dict[i_phase].density;

                                /* don't drift across interfaces; both cells must be cells the update acts on */
                                if((local_mass_c0 * local_mass_c1 > 0.0) &&
                                   (_C(i_layer).vof[i_frame_c0][c0][i_phase] > vof_cut) &&
                                   (_C(i_layer).vof[i_frame_c1][c1][i_phase] > vof_cut)){

                                    drift_volume = 0.0;

                                    for(i_dim = 0; i_dim < 3; i_dim++){

                                        /* drift_volume is defined per face, w/o information of vof's of c0/c1 */
                                        drift_volume +=   _F(i_layer).local_drift_adaption[i_face][i_data]
                                                        * rcfd->Data_Dict[i_phase][i_data].drift_velocity[i_dim]
                                                        * _F(i_layer).area[i_face][i_dim]
                                                        * rcfd->Solver->timestep_width[i_layer]
                                                        / (double)rcfd->Solver_Dict->number_of_drift_loops;
                                    }

                                    /* no flux through this face: nothing to exchange. Without this the
                                     * acceptor limiter below turned a zero flux into a negative exchange
                                     * whenever the cell taken as acceptor had a concentration above one,
                                     * and the two nodes of an ext face would pick different acceptors. */
                                    if(drift_volume == 0.0) continue;

                                    /* flip cells, such that flux is from c0 to c1 */
                                    if(drift_volume < 0.0){

                                        c = c0; c0 = c1; c1 = c;

                                        i_frame = i_frame_c0; i_frame_c0 = i_frame_c1; i_frame_c1 = i_frame;

                                        local_mass = local_mass_c0; local_mass_c0 = local_mass_c1; local_mass_c1 = local_mass;

                                        drift_volume *= -1.0;
                                    }

                                    /* in granular flow, drifting ability depend on solid volume fraction */
                                    if(rcfd->Phase_Dict[i_phase].hindered_drift_on == 1){

                                        hindering_factor = 1.0 - _C(i_layer).vof[i_frame_c0][c0][i_phase] / rcfd->Phase_Dict[i_phase].vof_max;

                                        if(hindering_factor < 0.0){

                                            hindering_factor = 0.0;
                                        }

                                        drift_volume *= hindering_factor;
                                    }
                                    /* TUHH 11/2026: no size drift in the dilute splash zone */
                                    else if((rcfd->Phase_Dict[i_phase].hindered_drift_on == 2) &&
                                            (rcfd->Phase_Dict[i_phase].drift_dilute_vof > 0.0)){

                                        hindering_factor = _C(i_layer).vof[i_frame_c0][c0][i_phase] / rcfd->Phase_Dict[i_phase].drift_dilute_vof;

                                        if(hindering_factor > 1.0){

                                            hindering_factor = 1.0;
                                        }

                                        drift_volume *= hindering_factor;
                                    }

                                    if(drift_volume > _C(i_layer).volume[c0]){

                                        drift_volume = _C(i_layer).volume[c0];
                                    }

                                    local_drift_exchange =    _C(i_layer).data[i_phase][c0][i_data]
                                                            * drift_volume
                                                            * _C(i_layer).vof[i_frame_c0][c0][i_phase]
                                                            * rcfd->Phase_Dict[i_phase].density;    /* [kg] */

                                    available_c1_mass   =     (1.0 - _C(i_layer).data[i_phase][c1][i_data])
                                                            * _C(i_layer).volume[c1]
                                                            * _C(i_layer).vof[i_frame_c1][c1][i_phase]
                                                            * rcfd->Phase_Dict[i_phase].density;

                                    if(local_drift_exchange > available_c1_mass){

                                        local_drift_exchange = available_c1_mass;
                                    }

                                    /* internal (owner) cells only, see the note at the top */
                                    if(c0 < int_cell_count) _C(i_layer).drift_exchange[c0] -= local_drift_exchange;
                                    if(c1 < int_cell_count) _C(i_layer).drift_exchange[c1] += local_drift_exchange;

                                    /* partition-boundary balance: what leaves the internal cell here must
                                     * enter the internal cell of the neighbour node; the global sum is zero */
                                    if(c1 >= int_cell_count){ seq_bnd_net += local_drift_exchange; seq_bnd_abs += 0.5 * fabs(local_drift_exchange); }
                                    if(c0 >= int_cell_count){ seq_bnd_net -= local_drift_exchange; seq_bnd_abs += 0.5 * fabs(local_drift_exchange); }
                                }
                            }

                            if(RCFD_DRIFT_DIAG_EVERY > 0 && (rcfd->Solver->global_run_counter % RCFD_DRIFT_DIAG_EVERY) == 0){

                                double g_bnd_net = PRF_GRSUM1(seq_bnd_net), g_bnd_abs = PRF_GRSUM1(seq_bnd_abs);

                                report_print0(rcfd, "[drift] step %d phase %d data %d loop %d/%d: partition boundaries exchange %.4e kg, imbalance %+.3e kg (must be round-off)",
                                              rcfd->Solver->global_run_counter, i_phase, i_data, i_drift + 1, (int)rcfd->Solver_Dict->number_of_drift_loops, g_bnd_abs, g_bnd_net);
                            }

                            /* no SUM_CELLS exchange (it double counted the partition-boundary fluxes) */

                            /* adapt _C(i_layer).data of internal cells by _C(i_layer).drift_exchange */
                            for(i_cell = 0; i_cell < int_cell_count; i_cell++){

                                local_mass =      _C(i_layer).volume[i_cell]
                                                * _C(i_layer).vof[i_frame][i_cell][i_phase]
                                                * rcfd->Phase_Dict[i_phase].density;

                                if(_C(i_layer).vof[i_frame][i_cell][i_phase] > rcfd->Solver_Dict->vof_cutoff_limit){

                                    _C(i_layer).data[i_phase][i_cell][i_data] = ((_C(i_layer).data[i_phase][i_cell][i_data] * local_mass) + _C(i_layer).drift_exchange[i_cell]) / local_mass;

                                    // adjust other conc's such that they still sum up to one
                                    if(rcfd->Phase_Dict[i_phase].drift_control_conc_sum_on)
                                    {
                                        current_data_index = i_data;

                                        current_conc = _C(i_layer).data[i_phase][i_cell][i_data];

                                        sum_of_other_concs = 0.0;

                                        for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

                                            if(i_data != current_data_index){

                                                sum_of_other_concs += _C(i_layer).data[i_phase][i_cell][i_data];
                                            }
                                        }

                                        if(sum_of_other_concs > 1.0e-6){

                                            conc_error = 1.0 - (sum_of_other_concs + current_conc);

                                            for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

                                                if(i_data != current_data_index){

                                                    _C(i_layer).data[i_phase][i_cell][i_data] += _C(i_layer).data[i_phase][i_cell][i_data] / sum_of_other_concs * conc_error;

                                                    if(_C(i_layer).data[i_phase][i_cell][i_data] < 0.0) _C(i_layer).data[i_phase][i_cell][i_data] = 0.0;
                                                    if(_C(i_layer).data[i_phase][i_cell][i_data] > 1.0) _C(i_layer).data[i_phase][i_cell][i_data] = 1.0;
                                                }
                                            }
                                        }

                                        i_data = current_data_index;
                                    }
                                }
                            }

                            /* corona copies follow their owner cells: this data and, through
                             * drift_control_conc_sum_on, the other concentrations changed */
                            if(rcfd->Solver_Dict->parallel_diffusion_on){

                                if(rcfd->Phase_Dict[i_phase].drift_control_conc_sum_on){

                                    parallel_map_phase_data(rcfd, i_phase);
                                }
                                else{

                                    parallel_communicate_data(rcfd, i_layer, i_phase, i_data, RCFD_PARALLEL_MAP_HOST_FORMAT);
                                }
                            }

                            break;
                        }
                        case RCFD_DATA_TEMPERATURE:
                        {
                            /* to do: replace local_drift_exchange by drift_exchange */
                            break;
                        }
                        default:
                            break;
                    }
                }
            }
        }
#undef LOCAL_DEBUG_FLAG

#endif
    }

    void swap_binarize_data(rCFD_type* rcfd, const int i_phase)
    {
#if RP_NODE

    /* Existing limitations
        L1: only single phase
    */

#define LOCAL_DEBUG_FLAG    0

#if 1   // local vars
        int     i_layer, i_data, i_cell, i_face;
        int     c0, c1, c;
        double  vol_flip;
#endif

        i_layer = rcfd->Solver->current_layer;

        swap_init(rcfd, i_phase);

        /* data as outer loop, because it is unlikely that many data get binarized */
        for(i_data = 0; i_data < rcfd->Phase_Dict[i_phase].number_of_data; i_data++){

            if(rcfd->Data_Dict[i_phase][i_data].type == RCFD_DATA_BINARY){

                /* Init data_swap */
                for(i_cell = 0; i_cell < _CELL_DICT(i_layer).number_of_cells; i_cell++){

                    _C(i_layer).data_swap[i_cell][i_data] = 0.0;
                }

                /* Artificial Diffusion (mass neutral) */
                for(i_face = 0; i_face < _FACE_DICT(i_layer).number_of_faces; i_face++){

                    c0 = _F(i_layer).c0[i_face];
                    c1 = _F(i_layer).c1[i_face];

                    if(_C(i_layer).data[i_phase][c0][i_data] > _C(i_layer).data[i_phase][c1][i_data]){

                        c = c0; c0 = c1; c1 = c;
                    }

                    if(_C(i_layer).data[i_phase][c1][i_data] > _C(i_layer).data[i_phase][c0][i_data]){

                        if(_C(i_layer).volume[c0] < _C(i_layer).volume[c1]){

                            vol_flip = _C(i_layer).volume[c0] * rcfd->Data_Dict[i_phase][i_data].binarization_art_diff;
                        }
                        else{

                            vol_flip = _C(i_layer).volume[c1] * rcfd->Data_Dict[i_phase][i_data].binarization_art_diff;
                        }

                        _C(i_layer).data_swap[c1][i_data] -=   vol_flip
                                                                      * (  _C(i_layer).data[i_phase][c1][i_data]
                                                                         - _C(i_layer).data[i_phase][c0][i_data])
                                                                      * 0.5;

                        _C(i_layer).data_swap[c0][i_data] +=   vol_flip
                                                                      * (  _C(i_layer).data[i_phase][c1][i_data]
                                                                         - _C(i_layer).data[i_phase][c0][i_data])
                                                                      * 0.5;
                    }
                }

                /* Binarization (mass acting) */
                for(i_cell = 0; i_cell < _CELL_DICT(i_layer).number_of_cells; i_cell++){

                    if(_C(i_layer).data[i_phase][i_cell][i_data] > 0.5){

                        _C(i_layer).data[i_phase][i_cell][i_data] = 1.0;
                    }
                    else{

                        _C(i_layer).data[i_phase][i_cell][i_data] = 0.0;
                    }
                }
            }
        }

#undef LOCAL_DEBUG_FLAG

#endif
    }

#endif
