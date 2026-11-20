/* rcfd_of_images.h - image cells behind periodic patches (TUHH 11/2026)
 *
 * Included by rcfd_of_bridge.c.
 *
 * The rCFD library has no neighbour across a periodic patch: it links a corona cell with the
 * cell it copies through the cell centre, and the copy of a cell on the other side of the
 * domain would need another centre. So the solver hands over, behind every face of a
 * periodic patch, an image cell: a copy of the cell on the other side, placed where that
 * cell would be if the mesh continued, and listed as a cell of this rank. The library then
 * interpolates the shifts across the periodic patch as inside the mesh.
 *
 * The library writes the shifts in its own cell numbering: cells of this rank, image cells,
 * corona cells. The run mode has no image cells. After the last frame, in every shift file
 * of this rank:
 *   target is an image cell      -> the shift is removed
 *   source is an image cell      -> the source becomes the cell the image copies: a cell of
 *                                   this rank, or a shift from the rank that holds it
 *   source is a corona cell      -> its index without the image cells
 *   incoming shift from an image cell of another rank -> from the cell that image copies
 *   incoming shift from a corona cell of another rank -> its index without the image cells
 * The vof files lose the lines of the image cells.
 */

static int      g_n_img = 0, g_n_fimg = 0, *g_img_host = NULL;

#define N_TRUE(T) ((T)->n_int - g_n_img)    /* cells of this rank without the image cells */

typedef struct { int c0, n0, c1; unsigned char w; } bridge_in_shift_t;

static void bridge_remove_images_from_database(int n_frames, int n_phases, int vof_files)
{
    long long img = g_n_img, g_img = 0;
    int     n_true = of_topo.n_int - g_n_img, n_lib = of_topo.n_int;
    int     i, p, r, i_frame;
    int     *g2l, *loc, *cnt, *ntrue, *disp, *own_host, *own_loc, *all_host, *all_loc, total = 0;
    long long removed = 0, to_own = 0, to_other = 0, sum[3], g_sum[3];
    const int wide = (of_size > 256);
    const size_t rec_a = 9, rec_b = wide ? 13 : 10;

    MPI_Allreduce(&img, &g_img, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if(g_img == 0) return;

    /* index of every cell on the rank that holds it */
    loc = (int*)calloc((size_t)of_n_cells_global, sizeof(int));
    g2l = (int*)calloc((size_t)of_n_cells_global, sizeof(int));

    if(loc == NULL || g2l == NULL) of_abort("out of memory (image cells)");

    for(i = 0; i < n_true; i++) loc[of_topo.gid[i]] = i;

    MPI_Allreduce(loc, g2l, of_n_cells_global, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    free(loc);

    /* image cells of all ranks: which cell they copy */
    cnt   = (int*)of_malloc((size_t)of_size * sizeof(int), "image counts");
    ntrue = (int*)of_malloc((size_t)of_size * sizeof(int), "cell counts");
    disp  = (int*)of_malloc((size_t)of_size * sizeof(int), "image offsets");

    MPI_Allgather(&g_n_img, 1, MPI_INT, cnt, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&n_true, 1, MPI_INT, ntrue, 1, MPI_INT, MPI_COMM_WORLD);

    for(r = 0; r < of_size; r++){ disp[r] = total; total += cnt[r]; }

    own_host = (int*)of_malloc((size_t)(g_n_img > 0 ? g_n_img : 1) * sizeof(int), "image owners");
    own_loc  = (int*)of_malloc((size_t)(g_n_img > 0 ? g_n_img : 1) * sizeof(int), "image owners");
    all_host = (int*)of_malloc((size_t)(total > 0 ? total : 1) * sizeof(int), "image owners");
    all_loc  = (int*)of_malloc((size_t)(total > 0 ? total : 1) * sizeof(int), "image owners");

    for(i = 0; i < g_n_img; i++){

        own_host[i] = g_img_host[i];
        own_loc[i]  = g2l[of_topo.gid[n_true + i]];
    }

    MPI_Allgatherv(own_host, g_n_img, MPI_INT, all_host, cnt, disp, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgatherv(own_loc,  g_n_img, MPI_INT, all_loc,  cnt, disp, MPI_INT, MPI_COMM_WORLD);

    free(g2l);

    for(p = 0; p < n_phases; p++){

        char    name[256], tmp[300];
        FILE    *f;
        long long size;
        unsigned char *buf, *q, *end;

        snprintf(name, sizeof(name), "%sdata/c2c/c2c_0_%d_%d", g_dir, p, of_rank);

        size = of_file_size(name);

        if(size < 0) of_abort("'%s' is missing", name);

        buf = (unsigned char*)of_malloc((size_t)(size > 0 ? size : 1), "shift file");

        f = fopen(name, "rb");

        if(f == NULL || fread(buf, 1, (size_t)size, f) != (size_t)size) of_abort("cannot read '%s'", name);

        fclose(f);

        snprintf(tmp, sizeof(tmp), "%s.tmp", name);

        f = fopen(tmp, "wb");

        if(f == NULL) of_abort("cannot write '%s': %s", tmp, strerror(errno));

        q = buf; end = buf + size;

        for(i_frame = 0; i_frame < n_frames; i_frame++){

            int     head[2], na, nb, k, n_a = 0, n_b = 0;
            unsigned char       *out_a;
            bridge_in_shift_t   *out_b;

            if(q + 8 > end) of_abort("%s: frame %d is missing", name, i_frame);

            memcpy(head, q, 8); q += 8;

            na = head[0]; nb = head[1];

            if(na < 0 || nb < 0 || q + (size_t)na * rec_a + (size_t)nb * rec_b > end) of_abort("%s: frame %d is damaged", name, i_frame);

            out_a = (unsigned char*)of_malloc(((size_t)na + (size_t)nb + 1) * rec_a, "shifts");
            out_b = (bridge_in_shift_t*)of_malloc(((size_t)na + (size_t)nb + 1) * sizeof(bridge_in_shift_t), "shifts");

            for(k = 0; k < na + nb; k++){

                int c0, c1, n0 = of_rank;
                unsigned char w;

                if(k < na)   { memcpy(&c0, q, 4); memcpy(&c1, q + 4, 4); w = q[8]; q += rec_a; }
                else if(wide){ memcpy(&c0, q, 4); memcpy(&n0, q + 4, 4); memcpy(&c1, q + 8, 4); w = q[12]; q += rec_b; }
                else         { memcpy(&c0, q, 4); n0 = (int)q[4]; memcpy(&c1, q + 5, 4); w = q[9]; q += rec_b; }

                /* target */
                if(c1 >= n_true && c1 < n_lib){ removed++; continue; }

                if(c1 < 0 || c1 >= n_lib) of_abort("%s: a shift ends in cell %d (cells of this rank: %d)", name, c1, n_true);

                /* source */
                if(k < na){

                    if(c0 >= n_lib) c0 -= g_n_img;                              /* corona cell */
                    else if(c0 >= n_true){                                      /* image cell of this rank */

                        n0 = own_host[c0 - n_true];
                        c0 = own_loc[c0 - n_true];

                        if(n0 == of_rank) to_own++; else to_other++;
                    }
                }
                else{

                    if(n0 < 0 || n0 >= of_size) of_abort("%s: incoming shift from rank %d", name, n0);

                    if(c0 < 0) of_abort("%s: incoming shift from cell %d of rank %d", name, c0, n0);

                    if(c0 >= ntrue[n0] + cnt[n0]) c0 -= cnt[n0];                /* corona cell of the other rank */
                    else if(c0 >= ntrue[n0]){                                   /* image cell of the other rank */

                        int j = disp[n0] + c0 - ntrue[n0];

                        n0 = all_host[j];
                        c0 = all_loc[j];

                        if(n0 == of_rank) to_own++; else to_other++;
                    }
                }

                if(n0 == of_rank){

                    unsigned char *o = out_a + (size_t)n_a * rec_a;

                    memcpy(o, &c0, 4); memcpy(o + 4, &c1, 4); o[8] = w;

                    n_a++;
                }
                else{

                    out_b[n_b].c0 = c0; out_b[n_b].n0 = n0; out_b[n_b].c1 = c1; out_b[n_b].w = w;

                    n_b++;
                }
            }

            /* incoming shifts grouped by rank, as the library writes them; the order within
               a rank is kept */
            {
                int sorted = 1;

                for(k = 1; k < n_b; k++){ if(out_b[k - 1].n0 > out_b[k].n0){ sorted = 0; break; } }

                if(!sorted){

                    bridge_in_shift_t *srt = (bridge_in_shift_t*)of_malloc(((size_t)n_b + 1) * sizeof(bridge_in_shift_t), "shifts");
                    int *start = (int*)calloc((size_t)of_size + 1, sizeof(int));

                    if(start == NULL) of_abort("out of memory (image cells)");

                    for(k = 0; k < n_b; k++) start[out_b[k].n0 + 1]++;
                    for(r = 0; r < of_size; r++) start[r + 1] += start[r];
                    for(k = 0; k < n_b; k++) srt[start[out_b[k].n0]++] = out_b[k];

                    memcpy(out_b, srt, (size_t)n_b * sizeof(bridge_in_shift_t));

                    free(srt); free(start);
                }
            }

            head[0] = n_a; head[1] = n_b;

            if(fwrite(head, 4, 2, f) != 2 || fwrite(out_a, rec_a, (size_t)n_a, f) != (size_t)n_a) of_abort("error while writing '%s'", tmp);

            for(k = 0; k < n_b; k++){

                unsigned char rec[13];

                memcpy(rec, &out_b[k].c0, 4);

                if(wide){ memcpy(rec + 4, &out_b[k].n0, 4); memcpy(rec + 8, &out_b[k].c1, 4); rec[12] = out_b[k].w; }
                else    { rec[4] = (unsigned char)out_b[k].n0; memcpy(rec + 5, &out_b[k].c1, 4); rec[9] = out_b[k].w; }

                if(fwrite(rec, rec_b, 1, f) != 1) of_abort("error while writing '%s'", tmp);
            }

            free(out_a); free(out_b);
        }

        if(q != end) of_abort("%s: %lld bytes after the last frame", name, (long long)(end - q));

        if(fclose(f) != 0) of_abort("error while writing '%s'", tmp);

        free(buf);

        if(rename(tmp, name) != 0) of_abort("cannot replace '%s': %s", name, strerror(errno));

        /* vof file (text): per frame the number of cells, then one value per cell of this rank */
        if(vof_files){

            FILE *fi, *fo;
            int  n_all, k;
            char line[128];

            snprintf(name, sizeof(name), "%sdata/vof/vof_0_%d_%d", g_dir, p, of_rank);
            snprintf(tmp, sizeof(tmp), "%s.tmp", name);

            fi = fopen(name, "r"); fo = fopen(tmp, "w");

            if(fi == NULL || fo == NULL) of_abort("cannot rewrite '%s'", name);

            for(i_frame = 0; i_frame < n_frames; i_frame++){

                if(fscanf(fi, "%d ", &n_all) != 1 || n_all != of_topo.n_cells) of_abort("%s: frame %d does not start with the number of cells", name, i_frame);

                fprintf(fo, "%d \n", n_all - g_n_img);

                for(k = 0; k < n_lib; k++){

                    if(fgets(line, sizeof(line), fi) == NULL) of_abort("%s: frame %d is too short", name, i_frame);

                    if(k < n_true) fputs(line, fo);
                }
            }

            fclose(fi);

            if(fclose(fo) != 0 || rename(tmp, name) != 0) of_abort("cannot replace '%s'", name);
        }
    }

    free(cnt); free(ntrue); free(disp); free(own_host); free(own_loc); free(all_host); free(all_loc);

    sum[0] = removed; sum[1] = to_own; sum[2] = to_other;

    MPI_Allreduce(sum, g_sum, 3, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    of_msg("image cells taken out of the database: %lld shifts into image cells removed; %lld shifts out of image cells now start at the cell behind the periodic patch (%lld of them on another rank)",
           g_sum[0], g_sum[1] + g_sum[2], g_sum[2]);
}

/* ------------------------------------------------------------------------- */
/* topology of the run mode                                                    */
/*                                                                             */
/* R = T without the image cells, plus the connections across the periodic     */
/* patches as ordinary faces:                                                  */
/*   both cells on this rank   -> a face between own cells (listed once)       */
/*   the other cell elsewhere  -> a corona cell (the copy of that cell, at its */
/*                                real position) and a partition face          */
/* The run mode then exchanges across periodic patches in the face swaps like  */
/* across any other face. The distance of the two cell centres across such a   */
/* face is not the distance of their positions: face index and distance go to  */
/* data/topo/periodic_<rank>.                                                  */
/* ------------------------------------------------------------------------- */

static int bridge_cmp_int_pair(const void *a, const void *b)
{
    const int *x = (const int*)a, *y = (const int*)b;

    return (x[0] > y[0]) - (x[0] < y[0]);
}

static void bridge_build_run_topology(const of_topo_t *T, of_topo_t *R, int **per_face, double **per_dist, int *n_per)
{
    const int n_true = N_TRUE(T), n_fint_true = T->n_fint - g_n_fimg;
    int     i, j, d, n_new = 0, n_fi = 0, n_fe = 0, n_cells, n_faces, k;
    int     *g2l = NULL, *img_c1, *img_kind, *cor, *new_gid, *new_host;
    double  *gx = NULL, *dist;
    long long img = g_n_img, g_img = 0;

    MPI_Allreduce(&img, &g_img, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    /* per face to an image cell */
    img_c1   = (int*)of_malloc((size_t)(g_n_fimg + 1) * sizeof(int), "periodic faces");
    img_kind = (int*)of_malloc((size_t)(g_n_fimg + 1) * sizeof(int), "periodic faces");     /* 0 dropped, 1 own cells, 2 partition */
    dist     = (double*)of_malloc((size_t)(g_n_fimg + 1) * sizeof(double), "periodic faces");
    new_gid  = (int*)of_malloc((size_t)(g_n_img + 1) * sizeof(int), "periodic corona");
    new_host = (int*)of_malloc((size_t)(g_n_img + 1) * sizeof(int), "periodic corona");

    if(g_img > 0){

        /* index on its rank, centre and volume of every cell of the mesh */
        size_t  n = (size_t)of_n_cells_global;
        int     *li = (int*)calloc(n, sizeof(int));
        double  *lx = (double*)calloc(4 * n, sizeof(double));

        g2l = (int*)calloc(n, sizeof(int));
        gx  = (double*)calloc(4 * n, sizeof(double));

        if(li == NULL || lx == NULL || g2l == NULL || gx == NULL) of_abort("out of memory (periodic faces)");

        for(i = 0; i < n_true; i++){

            li[T->gid[i]] = i;

            for(d = 0; d < 3; d++) lx[4 * (size_t)T->gid[i] + d] = T->x[3 * i + d];

            lx[4 * (size_t)T->gid[i] + 3] = T->vol[i];
        }

        MPI_Allreduce(li, g2l, of_n_cells_global, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(lx, gx, (int)(4 * n), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        free(li); free(lx);
    }

    /* corona cells by global index: (gid, index among the corona cells) */
    cor = (int*)of_malloc(2 * (size_t)(T->n_ext + g_n_img + 1) * sizeof(int), "corona index");

    for(i = 0; i < T->n_ext; i++){ cor[2 * i] = T->gid[T->n_int + i]; cor[2 * i + 1] = i; }

    qsort(cor, (size_t)T->n_ext, 2 * sizeof(int), bridge_cmp_int_pair);

    for(j = 0; j < g_n_fimg; j++){

        const int f  = n_fint_true + j;             /* face to an image cell */
        const int c0 = T->c0[f], ci = T->c1[f];
        int     gid, host;
        double  s2 = 0.0;

        if(ci < n_true || ci >= T->n_int || c0 < 0 || c0 >= n_true) of_abort("face %d (cells %d, %d) does not lead from a cell to an image cell", f, c0, ci);

        gid = T->gid[ci]; host = g_img_host[ci - n_true];

        for(d = 0; d < 3; d++) s2 += (T->x[3 * ci + d] - T->x[3 * c0 + d]) * (T->x[3 * ci + d] - T->x[3 * c0 + d]);

        dist[j] = sqrt(s2);

        if(host == of_rank){

            const int c1 = g2l[gid];

            /* both sides list the pair: keep it once */
            if(c0 < c1){ img_kind[j] = 1; img_c1[j] = c1; n_fi++; }
            else       { img_kind[j] = 0; img_c1[j] = -1; }
        }
        else{

            int key[2], *hit;

            key[0] = gid; key[1] = 0;

            hit = (int*)bsearch(key, cor, (size_t)T->n_ext, 2 * sizeof(int), bridge_cmp_int_pair);

            if(hit != NULL) img_c1[j] = n_true + hit[1];
            else{

                /* a new corona cell, unless an earlier periodic face made it already */
                for(k = 0; k < n_new; k++){ if(new_gid[k] == gid) break; }

                if(k == n_new){ new_gid[n_new] = gid; new_host[n_new] = host; n_new++; }

                img_c1[j] = n_true + T->n_ext + k;
            }

            img_kind[j] = 2; n_fe++;
        }
    }

    free(cor);

    n_cells = n_true + T->n_ext + n_new;
    n_faces = n_fint_true + n_fi + T->n_fext + n_fe;

    memset(R, 0, sizeof(*R));

    R->n_cells = n_cells; R->n_int = n_true; R->n_ext = T->n_ext + n_new;
    R->n_faces = n_faces; R->n_fint = n_fint_true + n_fi; R->n_fext = T->n_fext + n_fe;
    R->n_nodes = 0;

    R->x    = (double*)of_malloc(3 * (size_t)n_cells * sizeof(double), "cell centres");
    R->vol  = (double*)of_malloc((size_t)n_cells * sizeof(double), "cell volumes");
    R->host = (int*)   of_malloc((size_t)n_cells * sizeof(int), "cell hosts");
    R->gid  = (int*)   of_malloc((size_t)n_cells * sizeof(int), "global cell ids");
    R->c0   = (int*)   of_malloc((size_t)(n_faces + 1) * sizeof(int), "face c0");
    R->c1   = (int*)   of_malloc((size_t)(n_faces + 1) * sizeof(int), "face c1");
    R->area = (double*)of_malloc(3 * (size_t)(n_faces + 1) * sizeof(double), "face areas");

    *per_face = (int*)of_malloc((size_t)(n_fi + n_fe + 1) * sizeof(int), "periodic faces");
    *per_dist = (double*)of_malloc((size_t)(n_fi + n_fe + 1) * sizeof(double), "periodic faces");
    *n_per    = 0;

    /* cells: own, corona, new corona */
    for(i = 0, k = 0; i < T->n_cells; i++){

        if(i >= n_true && i < T->n_int) continue;

        for(d = 0; d < 3; d++) R->x[3 * k + d] = T->x[3 * i + d];

        R->vol[k] = T->vol[i]; R->host[k] = T->host[i]; R->gid[k] = T->gid[i];

        k++;
    }

    for(i = 0; i < n_new; i++, k++){

        for(d = 0; d < 3; d++) R->x[3 * k + d] = gx[4 * (size_t)new_gid[i] + d];

        R->vol[k] = gx[4 * (size_t)new_gid[i] + 3]; R->host[k] = new_host[i]; R->gid[k] = new_gid[i];
    }

    /* faces: between own cells, periodic between own cells, partition, periodic partition */
    k = 0;

    for(i = 0; i < n_fint_true; i++, k++){

        R->c0[k] = T->c0[i]; R->c1[k] = T->c1[i];

        for(d = 0; d < 3; d++) R->area[3 * k + d] = T->area[3 * i + d];
    }

    for(j = 0; j < g_n_fimg; j++){

        if(img_kind[j] != 1) continue;

        R->c0[k] = T->c0[n_fint_true + j]; R->c1[k] = img_c1[j];

        for(d = 0; d < 3; d++) R->area[3 * k + d] = T->area[3 * (n_fint_true + j) + d];

        (*per_face)[*n_per] = k; (*per_dist)[*n_per] = dist[j]; (*n_per)++;

        k++;
    }

    for(i = T->n_fint; i < T->n_faces; i++, k++){

        R->c0[k] = T->c0[i]; R->c1[k] = T->c1[i] - g_n_img;

        for(d = 0; d < 3; d++) R->area[3 * k + d] = T->area[3 * i + d];
    }

    for(j = 0; j < g_n_fimg; j++){

        if(img_kind[j] != 2) continue;

        R->c0[k] = T->c0[n_fint_true + j]; R->c1[k] = img_c1[j];

        for(d = 0; d < 3; d++) R->area[3 * k + d] = T->area[3 * (n_fint_true + j) + d];

        (*per_face)[*n_per] = k; (*per_dist)[*n_per] = dist[j]; (*n_per)++;

        k++;
    }

    if(k != n_faces) of_abort("run topology: %d faces written, %d expected", k, n_faces);

    free(img_c1); free(img_kind); free(dist); free(new_gid); free(new_host); free(g2l); free(gx);
}
