/*=============================================================================
 * rCFD Spray Distribution Module
 * --------------------------------
 * Gaussian radial + 1/z² axial liquid flux distribution for bottom-spray
 * two-fluid nozzles in fluidized bed granulation.
 *
 * Author:  Asif (TUHH) — generated with Claude assistance
 * Date:    April 2026
 * Context: SFB/CRC 1615 SMART Reactor — Sub-project C04
 *
 * USAGE:
 *   Replace the nozzle cone definitions, _growth_zone(), and the growth
 *   section in rCFD_user_jump_granulation.h with the code from this module.
 *   Integration points are clearly marked with [REPLACE] tags.
 *
 * PHYSICS:
 *   Liquid mass flux at position (r, z) in the spray cone:
 *
 *       phi(r,z) = (z_ref / z)^2 * exp( -r^2 / (2 * sigma(z)^2) )
 *
 *   where:
 *       r     = radial distance from cone axis
 *       z     = axial distance from nozzle tip (along cone axis)
 *       sigma = alpha(P_atm) * R(z)   [Gaussian spread width]
 *       R(z)  = r1 + (r2 - r1) * z/H  [cone radius at height z]
 *       alpha = alpha_0 / (1 + beta * P_atm)  [pressure coupling]
 *
 *   The weight w(r,z) is normalized globally so that the total deposited
 *   mass rate equals the specified spray rate * solid_content * efficiency.
 *=============================================================================*/

#ifndef RCFD_SPRAY_DISTRIBUTION_H
#define RCFD_SPRAY_DISTRIBUTION_H

/* =========================================================================
 * SECTION 1: SPRAY PROCESS PARAMETERS — SET THESE TO YOUR EXPERIMENT
 * ========================================================================= */

/* --- Spray rate --- */
#define _total_spray_rate_kg_per_hr     40.0        /* Total liquid spray rate [kg/hr] */
#define _number_of_nozzles              3
#define _spray_solid_content            0.30        /* Mass fraction of solid in spray liquid [-] */
#define _coating_efficiency             0.90        /* Fraction of solid that deposits on particles [-] */
#define _spray_weight_cap   1.0

/* Effective solid deposition rate per nozzle [kg/s]:
 * = (40/3600) * 0.30 * 0.80 / 3 = 0.000889 kg/s per nozzle
 * = total across 3 nozzles: 0.002667 kg/s
 * Compare: your old _global_growth_rate = 0.002 kg/s */
#define _effective_growth_rate  ((_total_spray_rate_kg_per_hr / 3600.0) * \
                                 _spray_solid_content * _coating_efficiency)

/* --- Nozzle geometry --- */
#define _spray_half_angle_deg           20.0        /* Spray cone half-angle [degrees] */
#define _atomization_pressure_bar       2.0         /* Atomization air pressure [bar] */
#define _nozzle_orifice_radius          0.005       /* Nozzle opening radius [m] (5 mm) */
#define _spray_penetration_height       0.1725      /* Axial spray height [m] */

/* --- Gaussian distribution parameters --- */
/* sigma(z) = alpha(P) * R(z)
 * alpha(P) = alpha_0 / (1 + beta * P_atm)
 *   P=1 bar: alpha = 0.55/(1+0.12*1) = 0.491  (broad, low pressure)
 *   P=2 bar: alpha = 0.55/(1+0.12*2) = 0.443  (moderate)
 *   P=3 bar: alpha = 0.55/(1+0.12*3) = 0.404  (tighter)
 *   P=5 bar: alpha = 0.55/(1+0.12*5) = 0.344  (very tight core)
 *
 * These values are empirical fits to two-fluid nozzle data from
 * Kramm et al. (2023) CIT and Spraying Systems Bulletin 459C. */
#define _alpha_0                        0.55
#define _alpha_beta                     0.12

/* Derived: sigma ratio at the configured pressure */
#define _sigma_ratio    (_alpha_0 / (1.0 + _alpha_beta * _atomization_pressure_bar))

/* Reference axial distance for 1/z^2 normalization [m].
 * Using a small offset from nozzle tip to avoid singularity at z=0. */
#define _z_ref          0.01

/* Minimum axial distance to avoid division by zero [m] */
#define _z_min_clip     0.005


/* =========================================================================
 * SECTION 2: AUTO-GENERATED NOZZLE CONE GEOMETRY
 * =========================================================================
 * [REPLACE] your old nozzle_cone_t definition and nozzle_cones[2] array */

#if RP_NODE

typedef struct {
    double point1[3];       /* Apex / nozzle tip position [m] */
    double point2[3];       /* Base center (top of spray cone) [m] */
    double radius1;         /* Radius at nozzle tip [m] */
    double radius2;         /* Radius at top of spray cone [m] */
    int    nozzle_id;       /* Nozzle identifier (1, 2, 3) */
    double cone_height;     /* Precomputed: ||point2 - point1|| [m] */
    double axis_unit[3];    /* Precomputed: unit vector along cone axis */
} spray_nozzle_t;

/* Build nozzle geometry from physical parameters.
 * Call once during initialization. */
static spray_nozzle_t spray_nozzles[_number_of_nozzles];

void _init_spray_nozzle_geometry(void)
{
    /* Nozzle tip positions (bottom spray, upward along Y) */
    double nozzle_tips[3][3] = {
        {0.0, 0.0300, 0.125},   /* Nozzle 1 */
        {0.0, 0.0300, 0.375},   /* Nozzle 2 */
        {0.0, 0.0300, 0.625}    /* Nozzle 3 */
    };

    double half_angle_rad = _spray_half_angle_deg * M_PI / 180.0;
    double H = _spray_penetration_height;
    double r2 = _nozzle_orifice_radius + H * tan(half_angle_rad);

    int i;
    for (i = 0; i < _number_of_nozzles; i++) {

        /* Tip position */
        spray_nozzles[i].point1[0] = nozzle_tips[i][0];
        spray_nozzles[i].point1[1] = nozzle_tips[i][1];
        spray_nozzles[i].point1[2] = nozzle_tips[i][2];

        /* Top of spray cone (tip + H in Y-direction) */
        spray_nozzles[i].point2[0] = nozzle_tips[i][0];
        spray_nozzles[i].point2[1] = nozzle_tips[i][1] + H;
        spray_nozzles[i].point2[2] = nozzle_tips[i][2];

        spray_nozzles[i].radius1     = _nozzle_orifice_radius;
        spray_nozzles[i].radius2     = r2;
        spray_nozzles[i].nozzle_id   = i + 1;
        spray_nozzles[i].cone_height = H;

        /* For bottom-spray upward: axis is purely +Y */
        spray_nozzles[i].axis_unit[0] = 0.0;
        spray_nozzles[i].axis_unit[1] = 1.0;
        spray_nozzles[i].axis_unit[2] = 0.0;
    }

    Message0("\n=== Auto-Generated Spray Nozzle Geometry ===\n");
    Message0("  Half-angle: %.1f deg | Pressure: %.1f bar | sigma_ratio: %.3f\n",
             _spray_half_angle_deg, _atomization_pressure_bar, _sigma_ratio);
    Message0("  Penetration height: %.4f m | Orifice radius: %.4f m\n", H, _nozzle_orifice_radius);
    Message0("  Cone top radius: %.4f m\n", r2);
    Message0("  Effective growth rate: %.6f kg/s (total across %d nozzles)\n",
             _effective_growth_rate, _number_of_nozzles);
    for (i = 0; i < _number_of_nozzles; i++) {
        Message0("  Nozzle %d: tip=(%.3f, %.3f, %.3f)  top=(%.3f, %.3f, %.3f)\n",
                 spray_nozzles[i].nozzle_id,
                 spray_nozzles[i].point1[0], spray_nozzles[i].point1[1], spray_nozzles[i].point1[2],
                 spray_nozzles[i].point2[0], spray_nozzles[i].point2[1], spray_nozzles[i].point2[2]);
    }
    Message0("=============================================\n");
}


/* =========================================================================
 * SECTION 3: SPRAY WEIGHT FUNCTION
 * =========================================================================
 * Returns the unnormalized Gaussian spray weight for a cell at position (x,y,z).
 * Returns 0.0 if the point is outside all nozzle cones.
 *
 * The weight combines:
 *   (a) 1/z^2 axial decay (spray thinning with distance)
 *   (b) Gaussian radial profile exp(-r^2 / (2*sigma^2))
 *   (c) sigma(z) = alpha(P) * R(z) — pressure-dependent concentration
 *
 * Also returns the nozzle_id via pointer (0 if outside). */

double _spray_weight(const double *pos, int *out_nozzle_id)
{
    int i;
    double best_weight = 0.0;
    int best_nozzle = 0;
 
    for (i = 0; i < _number_of_nozzles; i++) {
 
        const spray_nozzle_t *nz = &spray_nozzles[i];
 
        double dx = pos[0] - nz->point1[0];
        double dy = pos[1] - nz->point1[1];
        double dz = pos[2] - nz->point1[2];
 
        double z_proj = dx * nz->axis_unit[0] +
                        dy * nz->axis_unit[1] +
                        dz * nz->axis_unit[2];
 
        if (z_proj < 0.0 || z_proj > nz->cone_height) continue;
 
        double z_eff = (z_proj < _z_min_clip) ? _z_min_clip : z_proj;
 
        double t = z_proj / nz->cone_height;
        double R_at_z = nz->radius1 + t * (nz->radius2 - nz->radius1);
 
        double proj_x = nz->point1[0] + z_proj * nz->axis_unit[0];
        double proj_y = nz->point1[1] + z_proj * nz->axis_unit[1];
        double proj_z = nz->point1[2] + z_proj * nz->axis_unit[2];
 
        double rx = pos[0] - proj_x;
        double ry = pos[1] - proj_y;
        double rz = pos[2] - proj_z;
        double r_dist = sqrt(rx*rx + ry*ry + rz*rz);
 
        if (r_dist > 1.1 * R_at_z) continue;
 
        double sigma = _sigma_ratio * R_at_z;
        if (sigma < 1.0e-6) sigma = 1.0e-6;
 
        double gauss_radial = exp(-r_dist * r_dist / (2.0 * sigma * sigma));
        double axial_decay = (_z_ref * _z_ref) / (z_eff * z_eff);
 
        double w = axial_decay * gauss_radial;
 
        /* FIX: Cap weight to prevent outliers near nozzle tip */
        if (w > _spray_weight_cap) w = _spray_weight_cap;
 
        if (w > best_weight) {
            best_weight = w;
            best_nozzle = nz->nozzle_id;
        }
    }
 
    if (out_nozzle_id != NULL) {
        *out_nozzle_id = best_nozzle;
    }
 
    return best_weight;
}

/* Convenience: check if point is in any spray zone (replaces old _growth_zone) */
int _growth_zone_new(double *pos)
{
    int nid = 0;
    double w = _spray_weight(pos, &nid);
    return (w > 1.0e-10) ? 1 : 0;
}

/* Get nozzle ID (replaces old _get_nozzle_id) */
int _get_nozzle_id_new(double *pos)
{
    int nid = 0;
    _spray_weight(pos, &nid);
    return (nid > 0) ? nid : -1;
}

/* =========================================================================
 * SECTION 4: STATIC SPRAY WEIGHT CACHE
 * =========================================================================
 * Spray zone is fixed Eulerian geometry — compute once at init,
 * cache per rank, reuse at every timestep. Avoids rCFD's runtime
 * cell-center remapping corrupting geometric detection.
 * ========================================================================= */
static double  *_cached_spray_w    = NULL;   /* per-cell weight [int_cell_count] */
static int     *_cached_spray_nid  = NULL;   /* per-cell nozzle_id               */
static int      _cached_cell_count = 0;

void _build_spray_weight_cache(void *handle, int i_layer)
{
    int i_cell;
    int int_cell_count = rcfd_get_int_cell_count(handle, i_layer);
    double **x = rcfd_get_cell_center(handle, i_layer);

    /* Free old cache if any (shouldn't happen, but be safe) */
    if (_cached_spray_w)   { free(_cached_spray_w);   _cached_spray_w   = NULL; }
    if (_cached_spray_nid) { free(_cached_spray_nid); _cached_spray_nid = NULL; }

    _cached_spray_w   = (double*)calloc(int_cell_count, sizeof(double));
    _cached_spray_nid = (int*)   calloc(int_cell_count, sizeof(int));

    if (!_cached_spray_w || !_cached_spray_nid) {
        Message("ERROR [rank %d]: spray cache alloc failed (%d cells)\n",
                myid, int_cell_count);
        return;
    }
    _cached_cell_count = int_cell_count;

    int n_cone_cells = 0;
    double norm_check = 0.0;

    for (i_cell = 0; i_cell < int_cell_count; i_cell++) {
        int nid = 0;
        double w = _spray_weight(x[i_cell], &nid);
        _cached_spray_w[i_cell]   = w;
        _cached_spray_nid[i_cell] = nid;
        if (w > 1.0e-10) {
            n_cone_cells++;
            norm_check += w;
        }
    }

    Message("[RANK %d] Spray cache built: %d cells in cone, "
            "raw weight sum = %.4e\n",
            myid, n_cone_cells, norm_check);
}

/* =========================================================================
 * SECTION 7: DIAGNOSTIC — SPRAY WEIGHT MAP OUTPUT
 * =========================================================================
 * Optional: Call this once to write spray weights to a file for ParaView
 * visualization. Add to user_post or a one-time init block. */

void _write_spray_weight_map(void *handle)
{
    int i_layer = 0;
    int i_cell;
    int int_cell_count = rcfd_get_int_cell_count(handle, i_layer);
    double **x = rcfd_get_cell_center(handle, i_layer);
 
    /* Each node writes its OWN file (no race condition) */
    char fname[256];
    sprintf(fname, "./spray_weight_map_node%d.csv", myid);
 
    FILE *f = fopen(fname, "w");
    if (f == NULL) {
        Message("ERROR: node %d could not open %s\n", myid, fname);
        return;
    }
 
    /* Header on every node file (merge script will handle) */
    fprintf(f, "x,y,z,spray_weight,nozzle_id,z_proj_mm,r_mm,r_over_R\n");
 
    for (i_cell = 0; i_cell < int_cell_count; i_cell++) {
 
        int nid = 0;
        double w = _spray_weight(x[i_cell], &nid);
 
        if (w > 1.0e-10 && nid >= 1 && nid <= _number_of_nozzles) {
 
            /* Compute r and z for debugging columns */
            const spray_nozzle_t *nz = &spray_nozzles[nid - 1];
 
            double dx = x[i_cell][0] - nz->point1[0];
            double dy = x[i_cell][1] - nz->point1[1];
            double dz_c = x[i_cell][2] - nz->point1[2];
 
            double z_proj = dx * nz->axis_unit[0] +
                            dy * nz->axis_unit[1] +
                            dz_c * nz->axis_unit[2];
 
            double proj_x = nz->point1[0] + z_proj * nz->axis_unit[0];
            double proj_y = nz->point1[1] + z_proj * nz->axis_unit[1];
            double proj_z = nz->point1[2] + z_proj * nz->axis_unit[2];
 
            double rx = x[i_cell][0] - proj_x;
            double ry = x[i_cell][1] - proj_y;
            double rz = x[i_cell][2] - proj_z;
            double r_dist = sqrt(rx*rx + ry*ry + rz*rz);
 
            double t_p = z_proj / nz->cone_height;
            double R_at_z = nz->radius1 + t_p * (nz->radius2 - nz->radius1);
            double r_over_R = (R_at_z > 1.0e-6) ? r_dist / R_at_z : 0.0;
 
            /* Clean fprintf — NO trailing comma, fixed format */
            fprintf(f, "%.6f,%.6f,%.6f,%.6e,%d,%.3f,%.3f,%.4f\n",
                    x[i_cell][0],
                    x[i_cell][1],
                    x[i_cell][2],
                    w,
                    nid,
                    z_proj * 1000.0,
                    r_dist * 1000.0,
                    r_over_R);
        }
    }
 
    fclose(f);
 
    /* Node 0 writes the merge script */
    if (myid == 0) {
 
        FILE *fm = fopen("./merge_spray_map.sh", "w");
        if (fm != NULL) {
            fprintf(fm, "#!/bin/bash\n");
            fprintf(fm, "# Merge per-node spray weight CSVs\n");
            fprintf(fm, "head -1 spray_weight_map_node0.csv > spray_weight_map.csv\n");
            fprintf(fm, "for f in spray_weight_map_node*.csv; do\n");
            fprintf(fm, "  tail -n +2 \"$f\" >> spray_weight_map.csv\n");
            fprintf(fm, "done\n");
            fprintf(fm, "rm -f spray_weight_map_node*.csv\n");
            fprintf(fm, "wc -l spray_weight_map.csv\n");
            fprintf(fm, "echo 'Done. Load spray_weight_map.csv in ParaView.'\n");
            fclose(fm);
        }
 
        Message0("\n  Per-node CSVs written. Run: bash merge_spray_map.sh\n");
    }
}

#endif /* RP_NODE */

#endif /* RCFD_SPRAY_DISTRIBUTION_H */
