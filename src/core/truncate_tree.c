// Write a copy of the input VELOCIraptor augmented-tree halo catalog, keeping only
// halos whose FOF group Tvir >= 1e4K (the atomic-cooling threshold used by
// gas_cooling() in src/physics/cooling.c), with descendant links re-derived so that
// a dropped halo's lineage still resolves correctly to its next kept descendant,
// however many snapshots ahead that is.
//
// This runs in place of dracarys() (see meraxes.c) when Flag_WriteTruncatedTree is
// set. every snapshot is always loaded for this rank's assigned forests.before this.

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "meraxes.h"
#include "misc_tools.h"
#include "tree_flags.h"
#include "truncate_tree.h"
#include "virial_properties.h"

#define TVIR_CUT 1.0e4

typedef struct
{
  double mass_unit_to_internal; // (Mass_unit_to_solarmass / 1e10), matches read_halos-velociraptor.c
  int comoving_or_physical;
  int comoving_unit_flag;
  int cosmological_sim;
  double length_unit_to_kpc;
  double mass_unit_to_solarmass;
  char version[STRLEN];
  double velocity_unit_to_kms;
} tree_units_t;

static void read_source_units(tree_units_t* units)
{
  if (run_globals.mpi_rank == 0) {
    char fname[STRLEN * 2 + 8];
    switch (run_globals.params.TreesID) {
      case VELOCIRAPTOR_TREES:
        sprintf(fname, "%s/trees/%s", run_globals.params.SimulationDir, run_globals.params.CatalogFilePrefix);
        break;
      case VELOCIRAPTOR_TREES_AUG:
        sprintf(fname, "%s/augmented_trees/%s", run_globals.params.SimulationDir, run_globals.params.CatalogFilePrefix);
        break;
      default:
        mlog_error("write_truncated_tree: TreesID must be VELOCIRAPTOR_TREES or VELOCIRAPTOR_TREES_AUG (the "
                   "source catalog to truncate), not the truncated-tree ID itself.");
        ABORT(EXIT_FAILURE);
        break;
    }
    hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fd < 0) {
      mlog_error("write_truncated_tree: failed to open source file %s to read units.", fname);
      ABORT(EXIT_FAILURE);
    }
    H5LTget_attribute_int(fd, "Header/Units", "Comoving_or_Physical", &(units->comoving_or_physical));
    H5LTget_attribute_int(fd, "Header/Units", "Comoving_unit_flag", &(units->comoving_unit_flag));
    H5LTget_attribute_int(fd, "Header/Units", "Cosmological_Sim", &(units->cosmological_sim));
    H5LTget_attribute_double(fd, "Header/Units", "Length_unit_to_kpc", &(units->length_unit_to_kpc));
    H5LTget_attribute_double(fd, "Header/Units", "Mass_unit_to_solarmass", &(units->mass_unit_to_solarmass));
    H5LTget_attribute_double(fd, "Header/Units", "Velocity_unit_to_kms", &(units->velocity_unit_to_kms));
    memset(units->version, 0, STRLEN);
    H5LTget_attribute_string(fd, "Header/Units", "VERSION", units->version);
    H5Fclose(fd);
    units->mass_unit_to_internal = units->mass_unit_to_solarmass / 1.0e10;
  }
  MPI_Bcast(units, sizeof(tree_units_t), MPI_BYTE, 0, run_globals.mpi_comm);
}

// Returns true iff every halo in this FOF group should be kept (Tvir >= 1e4K,
// using the exact halo_type=1/AC formula that gas_cooling() itself uses).
static inline bool fof_group_is_kept(const fof_group_t* fof_group)
{
  double tvir = Vvir_to_Tvir((double)fof_group->Vvir, 1);
  return tvir >= TVIR_CUT;
}

int write_truncated_tree(void)
{
  // read_halos.c's preload checks (initialize_halo_storage() and read_halos()
  // itself) all treat Flag_WriteTruncatedTree the same as FlagInteractive/FlagMCMC,
  // so every snapshot is always preloaded (SnapshotHalo/SnapshotFOFGroup/
  // SnapshotTreesInfo) before we get here, regardless of what the param file's
  // own FlagInteractive/FlagMCMC values are.
  assert(run_globals.params.FlagInteractive || run_globals.params.FlagMCMC ||
         run_globals.params.Flag_WriteTruncatedTree);

  int last_snap = 0;
  for (int ii = 0; ii < run_globals.NOutputSnaps; ii++)
    if (run_globals.ListOutputSnaps[ii] > last_snap)
      last_snap = run_globals.ListOutputSnaps[ii];
  const int n_snaps = last_snap + 1;
  assert(run_globals.NStoreSnapshots == n_snaps);

  mlog("Writing Tvir>=1e4K-truncated tree...", MLOG_OPEN | MLOG_TIMERSTART);

  tree_units_t units;
  read_source_units(&units);
  const double hubble_h = run_globals.params.Hubble_h;

  // ---- Pass 1: mark keep/drop per halo, per snapshot, and assign compact
  // ---- (per-rank-local) new indices to kept halos.
  bool** keep = malloc(sizeof(bool*) * n_snaps);
  int** new_index = malloc(sizeof(int*) * n_snaps);
  int* local_n_halos_kept = calloc(n_snaps, sizeof(int));
  int* local_n_fof_kept = calloc(n_snaps, sizeof(int));

  for (int s = 0; s < n_snaps; s++) {
    const int n_halos = run_globals.SnapshotTreesInfo[s].n_halos;
    const int n_fof = run_globals.SnapshotTreesInfo[s].n_fof_groups;
    halo_t* halos = run_globals.SnapshotHalo[s];
    fof_group_t* fof_groups = run_globals.SnapshotFOFGroup[s];

    keep[s] = (n_halos > 0) ? malloc(sizeof(bool) * n_halos) : NULL;
    new_index[s] = (n_halos > 0) ? malloc(sizeof(int) * n_halos) : NULL;
    if (n_halos <= 0)
      continue;

    bool* group_kept = malloc(sizeof(bool) * (n_fof > 0 ? n_fof : 1));
    for (int j = 0; j < n_fof; j++)
      group_kept[j] = fof_group_is_kept(&fof_groups[j]);

    int running_halos = 0;
    int running_fof = 0;
    for (int i = 0; i < n_halos; i++) {
      halo_t* h = &halos[i];
      int fof_idx = (int)(h->FOFGroup - fof_groups);
      if (fof_idx < 0 || fof_idx >= n_fof) {
        mlog_error("write_truncated_tree: snapshot %d halo %d has an out-of-range FOFGroup pointer "
                   "(fof_idx=%d, n_fof=%d) -- halo->FOFGroup does not point into this snapshot's FOF array.",
                   s,
                   i,
                   fof_idx,
                   n_fof);
        ABORT(EXIT_FAILURE);
      }
      bool kept = group_kept[fof_idx];
      keep[s][i] = kept;
      if (kept) {
        new_index[s][i] = running_halos++;
        if (h->Type == 0)
          running_fof++;
      } else {
        new_index[s][i] = -1;
      }
    }
    local_n_halos_kept[s] = running_halos;
    local_n_fof_kept[s] = running_fof;
    free(group_kept);
  }

  // ---- Pass 2: for every kept halo, walk DescIndex/SnapOffset forward (this
  // ---- rank's own in-memory data only -- forests never cross ranks) until we
  // ---- find another kept halo, or run out of tree. Store the *local* target
  // ---- snapshot/index (converted to a global new ID once we know per-snapshot
  // ---- MPI offsets, below).
  int** head_target_snap = malloc(sizeof(int*) * n_snaps);
  int** head_target_index = malloc(sizeof(int*) * n_snaps);

  for (int s = 0; s < n_snaps; s++) {
    const int n_halos = run_globals.SnapshotTreesInfo[s].n_halos;
    head_target_snap[s] = (n_halos > 0) ? malloc(sizeof(int) * n_halos) : NULL;
    head_target_index[s] = (n_halos > 0) ? malloc(sizeof(int) * n_halos) : NULL;

    for (int i = 0; i < n_halos; i++) {
      if (!keep[s][i])
        continue;

      int cur_s = s;
      int cur_i = i;
      int target_s = -1;
      int target_i = -1;
      // A halo's descendant chain strictly advances in snapshot number, so this
      // can loop at most n_snaps times; the guard is only a defensive backstop
      // against malformed tree data.
      for (int guard = 0; guard <= n_snaps + 1; guard++) {
        halo_t* h = &run_globals.SnapshotHalo[cur_s][cur_i];
        if (h->DescIndex == -1)
          break; // life ends here -- no kept descendant
        int next_s = cur_s + h->SnapOffset;
        if (next_s >= n_snaps || next_s <= cur_s)
          break; // out of range, or a non-advancing offset -- treat as terminal

        // h->DescIndex is the ORIGINAL (file-space) local index at next_s -- with
        // multiple ranks each holding only its own forests' halos, that has to be
        // translated to this rank's compacted runtime array position first, exactly
        // like the normal galaxy/halo linking code in dracarys.c does.
        int next_i = h->DescIndex;
        int* lookup = run_globals.SnapshotIndexLookup[next_s];
        if (lookup != NULL) {
          next_i = find_original_index(h->DescIndex, lookup, run_globals.SnapshotTreesInfo[next_s].n_halos);
          if (next_i == -1)
            break; // descendant isn't part of this rank's (this forest's) kept set -- treat as terminal
        }
        if (next_i < 0 || next_i >= run_globals.SnapshotTreesInfo[next_s].n_halos) {
          mlog_error("write_truncated_tree: snapshot %d halo %d's descendant chain points to an "
                     "out-of-range halo at snapshot %d (DescIndex=%d, translated index=%d, n_halos=%d).",
                     s,
                     i,
                     next_s,
                     h->DescIndex,
                     next_i,
                     run_globals.SnapshotTreesInfo[next_s].n_halos);
          ABORT(EXIT_FAILURE);
        }
        if (keep[next_s][next_i]) {
          target_s = next_s;
          target_i = new_index[next_s][next_i];
          break;
        }
        cur_s = next_s;
        cur_i = next_i;
      }
      head_target_snap[s][i] = target_s;
      head_target_index[s][i] = target_i;
    }
  }

  // ---- Determine per-snapshot, per-rank write offsets and global totals.
  // (Only halo-level offsets are needed for writing -- FOF-group membership in
  // the output schema is encoded via hostHaloID, not a separate indexed array.)
  int* halo_offset = malloc(sizeof(int) * n_snaps);
  int* global_n_halos = malloc(sizeof(int) * n_snaps);
  int* global_n_fof = malloc(sizeof(int) * n_snaps);
  int global_n_halos_max = 0;
  int global_n_fof_max = 0;

  for (int s = 0; s < n_snaps; s++) {
    int excl_halo = 0;
    MPI_Exscan(&local_n_halos_kept[s], &excl_halo, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    if (run_globals.mpi_rank == 0)
      excl_halo = 0;
    halo_offset[s] = excl_halo;

    MPI_Allreduce(&local_n_halos_kept[s], &global_n_halos[s], 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(&local_n_fof_kept[s], &global_n_fof[s], 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    if (global_n_halos[s] > global_n_halos_max)
      global_n_halos_max = global_n_halos[s];
    if (global_n_fof[s] > global_n_fof_max)
      global_n_fof_max = global_n_fof[s];

    // The *original* (pre-cut) counts also need to be summed across ranks --
    // SnapshotTreesInfo[s] only holds this rank's own share of the source catalog.
    int orig_n_halos = run_globals.SnapshotTreesInfo[s].n_halos;
    int orig_n_fof = run_globals.SnapshotTreesInfo[s].n_fof_groups;
    int global_orig_n_halos, global_orig_n_fof;
    MPI_Allreduce(&orig_n_halos, &global_orig_n_halos, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);
    MPI_Allreduce(&orig_n_fof, &global_orig_n_fof, 1, MPI_INT, MPI_SUM, run_globals.mpi_comm);

    mlog("snapshot %d :: kept %d/%d halos, %d/%d FOF groups",
         MLOG_MESG,
         s,
         global_n_halos[s],
         global_orig_n_halos,
         global_n_fof[s],
         global_orig_n_fof);
  }

  // ---- Create the output tree file + per-snapshot groups/datasets (rank 0 only).
  char tree_dir[STRLEN + 20];
  sprintf(tree_dir, "%s/truncated_trees", run_globals.params.SimulationDir);
  char tree_fname[STRLEN * 2 + 8];
  sprintf(tree_fname, "%s/%s", tree_dir, run_globals.params.CatalogFilePrefix);
  char stats_fname[STRLEN + 40];
  sprintf(stats_fname, "%s/meraxes_augmented_stats.h5", tree_dir);

  const char* dset_names_float[] = { "Mass_200crit", "Mass_tot", "R_200crit", "Vmax",
                                      "Xc",           "Yc",       "Zc",        "VXc",
                                      "VYc",          "VZc",      "AngMom" };
  const int n_dset_float = sizeof(dset_names_float) / sizeof(dset_names_float[0]);

  if (run_globals.mpi_rank == 0) {
    struct stat st;
    if (stat(tree_dir, &st) != 0)
      mkdir(tree_dir, 02755);

    hid_t fd = H5Fcreate(tree_fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (fd < 0) {
      mlog_error("write_truncated_tree: failed to create %s", tree_fname);
      ABORT(EXIT_FAILURE);
    }

    hid_t header_grp = H5Gcreate(fd, "Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t units_grp = H5Gcreate(header_grp, "Units", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5LTset_attribute_int(header_grp, "Units", "Comoving_or_Physical", &units.comoving_or_physical, 1);
    H5LTset_attribute_int(header_grp, "Units", "Comoving_unit_flag", &units.comoving_unit_flag, 1);
    H5LTset_attribute_int(header_grp, "Units", "Cosmological_Sim", &units.cosmological_sim, 1);
    H5LTset_attribute_double(header_grp, "Units", "Length_unit_to_kpc", &units.length_unit_to_kpc, 1);
    H5LTset_attribute_double(header_grp, "Units", "Mass_unit_to_solarmass", &units.mass_unit_to_solarmass, 1);
    H5LTset_attribute_string(header_grp, "Units", "VERSION", units.version);
    H5LTset_attribute_double(header_grp, "Units", "Velocity_unit_to_kms", &units.velocity_unit_to_kms, 1);
    H5Gclose(units_grp);
    H5Gclose(header_grp);

    for (int s = 0; s < n_snaps; s++) {
      char grp_name[16];
      sprintf(grp_name, "Snap_%03d", s);
      hid_t grp = H5Gcreate(fd, grp_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      int nhalos_s = global_n_halos[s];
      double scale_factor = 1.0 / (1.0 + run_globals.ZZ[s]);
      H5LTset_attribute_int(fd, grp_name, "NHalos", &nhalos_s, 1);
      H5LTset_attribute_int(fd, grp_name, "Snapnum", &s, 1);
      H5LTset_attribute_double(fd, grp_name, "scalefactor", &scale_factor, 1);

      hsize_t dim[1] = { (hsize_t)nhalos_s };
      hid_t dspace = H5Screate_simple(1, dim, NULL);

      hid_t dset;
      dset = H5Dcreate(grp, "ID", H5T_NATIVE_LONG, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose(dset);
      dset = H5Dcreate(grp, "Head", H5T_NATIVE_LONG, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose(dset);
      dset = H5Dcreate(grp, "hostHaloID", H5T_NATIVE_LONG, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose(dset);
      dset = H5Dcreate(grp, "ForestID", H5T_NATIVE_ULONG, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose(dset);
      dset = H5Dcreate(grp, "npart", H5T_NATIVE_UINT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose(dset);
      for (int k = 0; k < n_dset_float; k++) {
        dset = H5Dcreate(grp, dset_names_float[k], H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset);
      }

      H5Sclose(dspace);
      H5Gclose(grp);
    }
    H5Fclose(fd);
  }
  MPI_Barrier(run_globals.mpi_comm);

  // ---- Each rank writes its own disjoint slice of kept halos, per snapshot.
  hid_t plist_acc = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_acc, run_globals.mpi_comm, MPI_INFO_NULL);
  hid_t fd = H5Fopen(tree_fname, H5F_ACC_RDWR, plist_acc);
  H5Pclose(plist_acc);
  if (fd < 0) {
    mlog_error("write_truncated_tree: rank %d failed to reopen %s for writing", run_globals.mpi_rank, tree_fname);
    ABORT(EXIT_FAILURE);
  }

  hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_INDEPENDENT);

  for (int s = 0; s < n_snaps; s++) {
    const int n_local = local_n_halos_kept[s];
    char grp_name[16];
    sprintf(grp_name, "Snap_%03d", s);
    hid_t grp = H5Gopen(fd, grp_name, H5P_DEFAULT);

    if (n_local > 0) {
      int64_t* out_ID = malloc(sizeof(int64_t) * n_local);
      int64_t* out_Head = malloc(sizeof(int64_t) * n_local);
      int64_t* out_hostHaloID = malloc(sizeof(int64_t) * n_local);
      uint64_t* out_ForestID = malloc(sizeof(uint64_t) * n_local);
      uint32_t* out_npart = malloc(sizeof(uint32_t) * n_local);
      float* out_float[11];
      for (int k = 0; k < n_dset_float; k++)
        out_float[k] = malloc(sizeof(float) * n_local);

      const int n_halos = run_globals.SnapshotTreesInfo[s].n_halos;
      halo_t* halos = run_globals.SnapshotHalo[s];
      const double scale_factor = 1.0 / (1.0 + run_globals.ZZ[s]);
      int64_t last_host_new_id = -1;
      int row = 0;

      for (int i = 0; i < n_halos; i++) {
        if (!keep[s][i])
          continue;

        halo_t* h = &halos[i];
        int64_t self_id = (int64_t)s * 1000000000000LL + (int64_t)(halo_offset[s] + new_index[s][i]) + 1LL;

        out_ID[row] = self_id;
        out_ForestID[row] = (uint64_t)h->ForestID;
        out_npart[row] = (uint32_t)h->Len;

        int target_s = head_target_snap[s][i];
        if (target_s >= 0) {
          int64_t target_id =
            (int64_t)target_s * 1000000000000LL + (int64_t)(halo_offset[target_s] + head_target_index[s][i]) + 1LL;
          out_Head[row] = target_id;
        } else {
          out_Head[row] = self_id; // matches the reader's cyclic "life ends here" convention
        }

        if (h->Type == 0) {
          out_hostHaloID[row] = -1;
          last_host_new_id = self_id;

          bool below_thresh = (h->TreeFlags & TREE_CASE_BELOW_VIRIAL_THRESHOLD) != 0;
          fof_group_t* fof_group = h->FOFGroup;
          if (below_thresh) {
            out_float[0][row] = 0.0f;  // Mass_200crit
            out_float[2][row] = -1.0f; // R_200crit
          } else {
            out_float[0][row] =
              (float)((double)(fof_group->Mvir) / (double)(fof_group->FOFMvirModifier) /
                      (hubble_h * units.mass_unit_to_internal));
            out_float[2][row] = (float)((double)(fof_group->Rvir) / hubble_h);
          }
        } else {
          out_hostHaloID[row] = last_host_new_id;
          out_float[0][row] = (float)((double)(h->Mvir) / (hubble_h * units.mass_unit_to_internal)); // placeholder
          out_float[2][row] = 0.0f;                                                                   // placeholder
        }

        out_float[1][row] = (float)((double)(h->Mvir) / (hubble_h * units.mass_unit_to_internal)); // Mass_tot
        out_float[3][row] = h->Vmax;
        out_float[4][row] = (float)((double)(h->Pos[0]) * scale_factor / hubble_h); // Xc
        out_float[5][row] = (float)((double)(h->Pos[1]) * scale_factor / hubble_h); // Yc
        out_float[6][row] = (float)((double)(h->Pos[2]) * scale_factor / hubble_h); // Zc
        out_float[7][row] = (float)((double)(h->Vel[0]) * scale_factor);            // VXc
        out_float[8][row] = (float)((double)(h->Vel[1]) * scale_factor);            // VYc
        out_float[9][row] = (float)((double)(h->Vel[2]) * scale_factor);            // VZc
        out_float[10][row] = (float)((double)(h->AngMom) / hubble_h);               // AngMom

        row++;
      }
      assert(row == n_local);

      hsize_t offset[1] = { (hsize_t)halo_offset[s] };
      hsize_t count[1] = { (hsize_t)n_local };
      hid_t memspace = H5Screate_simple(1, count, NULL);

#define WRITE_COL(name, h5type, buf)                                                                                 \
  {                                                                                                                   \
    hid_t dset_id = H5Dopen(grp, name, H5P_DEFAULT);                                                                  \
    hid_t fspace_id = H5Dget_space(dset_id);                                                                          \
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);                                        \
    H5Dwrite(dset_id, h5type, memspace, fspace_id, xfer_plist, buf);                                                  \
    H5Sclose(fspace_id);                                                                                              \
    H5Dclose(dset_id);                                                                                                \
  }

      WRITE_COL("ID", H5T_NATIVE_LONG, out_ID);
      WRITE_COL("Head", H5T_NATIVE_LONG, out_Head);
      WRITE_COL("hostHaloID", H5T_NATIVE_LONG, out_hostHaloID);
      WRITE_COL("ForestID", H5T_NATIVE_ULONG, out_ForestID);
      WRITE_COL("npart", H5T_NATIVE_UINT, out_npart);
      for (int k = 0; k < n_dset_float; k++)
        WRITE_COL(dset_names_float[k], H5T_NATIVE_FLOAT, out_float[k]);

#undef WRITE_COL

      H5Sclose(memspace);
      free(out_ID);
      free(out_Head);
      free(out_hostHaloID);
      free(out_ForestID);
      free(out_npart);
      for (int k = 0; k < n_dset_float; k++)
        free(out_float[k]);
    }

    H5Gclose(grp);
  }

  H5Pclose(xfer_plist);
  H5Fclose(fd);
  MPI_Barrier(run_globals.mpi_comm);

  // ---- Regenerate the forests stats file. This rank's owned forest list is
  // ---- whatever it was assigned to read (run_globals.RequestedForestId), or --
  // ---- for a run with no forest partitioning at all (e.g. a single-rank test
  // ---- run) -- the distinct ForestIDs actually seen among this rank's kept
  // ---- halos.
  long* my_forest_ids = NULL;
  int n_my_forests = 0;
  bool owns_forest_id_list = false;

  if (run_globals.RequestedForestId != NULL && run_globals.NRequestedForests > 0) {
    n_my_forests = run_globals.NRequestedForests;
    my_forest_ids = malloc(sizeof(long) * n_my_forests);
    memcpy(my_forest_ids, run_globals.RequestedForestId, sizeof(long) * n_my_forests);
  } else {
    int cap = 1024;
    int n_seen = 0;
    long* seen = malloc(sizeof(long) * cap);
    for (int s = 0; s < n_snaps; s++) {
      const int n_halos = run_globals.SnapshotTreesInfo[s].n_halos;
      halo_t* halos = run_globals.SnapshotHalo[s];
      for (int i = 0; i < n_halos; i++) {
        if (!keep[s][i])
          continue;
        if (n_seen == cap) {
          cap *= 2;
          seen = realloc(seen, sizeof(long) * cap);
        }
        seen[n_seen++] = (long)halos[i].ForestID;
      }
    }
    qsort(seen, n_seen, sizeof(long), compare_longs);
    int n_unique = 0;
    for (int ii = 0; ii < n_seen; ii++)
      if (ii == 0 || seen[ii] != seen[ii - 1])
        seen[n_unique++] = seen[ii];
    my_forest_ids = seen;
    n_my_forests = n_unique;
    owns_forest_id_list = true;
  }
  (void)owns_forest_id_list;

  int* forest_n_halos = calloc(n_my_forests > 0 ? n_my_forests : 1, sizeof(int));
  int* forest_n_fof = calloc(n_my_forests > 0 ? n_my_forests : 1, sizeof(int));
  int* forest_max_contemp_halos = calloc(n_my_forests > 0 ? n_my_forests : 1, sizeof(int));
  int* forest_max_contemp_fof = calloc(n_my_forests > 0 ? n_my_forests : 1, sizeof(int));
  int** forest_snap_halos = malloc(sizeof(int*) * n_snaps);

  for (int s = 0; s < n_snaps; s++) {
    forest_snap_halos[s] = calloc(n_my_forests > 0 ? n_my_forests : 1, sizeof(int));
    int* snap_fof_count = calloc(n_my_forests > 0 ? n_my_forests : 1, sizeof(int));
    const int n_halos = run_globals.SnapshotTreesInfo[s].n_halos;
    halo_t* halos = run_globals.SnapshotHalo[s];

    for (int i = 0; i < n_halos; i++) {
      if (!keep[s][i])
        continue;
      long fid = (long)halos[i].ForestID;
      long* pos = bsearch(&fid, my_forest_ids, (size_t)n_my_forests, sizeof(long), compare_longs);
      if (pos == NULL)
        continue; // should not happen
      int fp = (int)(pos - my_forest_ids);
      forest_snap_halos[s][fp]++;
      forest_n_halos[fp]++;
      if (halos[i].Type == 0) {
        snap_fof_count[fp]++;
        forest_n_fof[fp]++;
      }
    }
    for (int fp = 0; fp < n_my_forests; fp++) {
      if (forest_snap_halos[s][fp] > forest_max_contemp_halos[fp])
        forest_max_contemp_halos[fp] = forest_snap_halos[s][fp];
      if (snap_fof_count[fp] > forest_max_contemp_fof[fp])
        forest_max_contemp_fof[fp] = snap_fof_count[fp];
    }
    free(snap_fof_count);
  }

  int* recvcounts = malloc(sizeof(int) * run_globals.mpi_size);
  MPI_Gather(&n_my_forests, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, run_globals.mpi_comm);
  int* displs = malloc(sizeof(int) * run_globals.mpi_size);
  int n_forests_global = 0;
  if (run_globals.mpi_rank == 0) {
    displs[0] = 0;
    for (int r = 0; r < run_globals.mpi_size; r++) {
      if (r > 0)
        displs[r] = displs[r - 1] + recvcounts[r - 1];
      n_forests_global += recvcounts[r];
    }
  }
  MPI_Bcast(&n_forests_global, 1, MPI_INT, 0, run_globals.mpi_comm);

  long* global_forest_ids = (run_globals.mpi_rank == 0) ? malloc(sizeof(long) * n_forests_global) : NULL;
  int* global_forest_n_halos = (run_globals.mpi_rank == 0) ? malloc(sizeof(int) * n_forests_global) : NULL;
  int* global_forest_n_fof = (run_globals.mpi_rank == 0) ? malloc(sizeof(int) * n_forests_global) : NULL;
  int* global_forest_max_halos = (run_globals.mpi_rank == 0) ? malloc(sizeof(int) * n_forests_global) : NULL;
  int* global_forest_max_fof = (run_globals.mpi_rank == 0) ? malloc(sizeof(int) * n_forests_global) : NULL;

  MPI_Gatherv(
    my_forest_ids, n_my_forests, MPI_LONG, global_forest_ids, recvcounts, displs, MPI_LONG, 0, run_globals.mpi_comm);
  MPI_Gatherv(forest_n_halos,
              n_my_forests,
              MPI_INT,
              global_forest_n_halos,
              recvcounts,
              displs,
              MPI_INT,
              0,
              run_globals.mpi_comm);
  MPI_Gatherv(
    forest_n_fof, n_my_forests, MPI_INT, global_forest_n_fof, recvcounts, displs, MPI_INT, 0, run_globals.mpi_comm);
  MPI_Gatherv(forest_max_contemp_halos,
              n_my_forests,
              MPI_INT,
              global_forest_max_halos,
              recvcounts,
              displs,
              MPI_INT,
              0,
              run_globals.mpi_comm);
  MPI_Gatherv(forest_max_contemp_fof,
              n_my_forests,
              MPI_INT,
              global_forest_max_fof,
              recvcounts,
              displs,
              MPI_INT,
              0,
              run_globals.mpi_comm);

  int** global_forest_snap_halos = NULL;
  if (run_globals.mpi_rank == 0) {
    global_forest_snap_halos = malloc(sizeof(int*) * n_snaps);
    for (int s = 0; s < n_snaps; s++)
      global_forest_snap_halos[s] = malloc(sizeof(int) * n_forests_global);
  }
  for (int s = 0; s < n_snaps; s++)
    MPI_Gatherv(forest_snap_halos[s],
                n_my_forests,
                MPI_INT,
                run_globals.mpi_rank == 0 ? global_forest_snap_halos[s] : NULL,
                recvcounts,
                displs,
                MPI_INT,
                0,
                run_globals.mpi_comm);

  if (run_globals.mpi_rank == 0) {
    hid_t sfd = H5Fcreate(stats_fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (sfd < 0) {
      mlog_error("write_truncated_tree: failed to create %s", stats_fname);
      ABORT(EXIT_FAILURE);
    }
    H5LTset_attribute_int(sfd, "/", "n_snaps", &n_snaps, 1);
    H5LTset_attribute_int(sfd, "/", "n_halos_max", &global_n_halos_max, 1);
    H5LTset_attribute_int(sfd, "/", "n_fof_groups_max", &global_n_fof_max, 1);

    hsize_t dim_snap[1] = { (hsize_t)n_snaps };
    H5LTmake_dataset(sfd, "n_halos", 1, dim_snap, H5T_NATIVE_INT, global_n_halos);
    H5LTmake_dataset(sfd, "n_fof_groups", 1, dim_snap, H5T_NATIVE_INT, global_n_fof);

    hid_t forests_grp = H5Gcreate(sfd, "forests", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5LTset_attribute_int(sfd, "forests", "n_forests", &n_forests_global, 1);
    hsize_t dim_forest[1] = { (hsize_t)n_forests_global };
    H5LTmake_dataset(forests_grp, "forest_ids", 1, dim_forest, H5T_NATIVE_LONG, global_forest_ids);
    H5LTmake_dataset(forests_grp, "n_halos", 1, dim_forest, H5T_NATIVE_INT, global_forest_n_halos);
    H5LTmake_dataset(forests_grp, "n_fof_groups", 1, dim_forest, H5T_NATIVE_INT, global_forest_n_fof);
    H5LTmake_dataset(forests_grp, "max_contemporaneous_halos", 1, dim_forest, H5T_NATIVE_INT, global_forest_max_halos);
    H5LTmake_dataset(
      forests_grp, "max_contemporaneous_fof_groups", 1, dim_forest, H5T_NATIVE_INT, global_forest_max_fof);
    H5Gclose(forests_grp);

    hid_t snaps_grp = H5Gcreate(sfd, "snapshots", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (int s = 0; s < n_snaps; s++) {
      char dset_name[16];
      sprintf(dset_name, "Snap%03d", s);
      H5LTmake_dataset(snaps_grp, dset_name, 1, dim_forest, H5T_NATIVE_INT, global_forest_snap_halos[s]);
    }
    H5Gclose(snaps_grp);
    H5Fclose(sfd);

    free(global_forest_ids);
    free(global_forest_n_halos);
    free(global_forest_n_fof);
    free(global_forest_max_halos);
    free(global_forest_max_fof);
    for (int s = 0; s < n_snaps; s++)
      free(global_forest_snap_halos[s]);
    free(global_forest_snap_halos);
  }

  MPI_Barrier(run_globals.mpi_comm);
  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

  // ---- cleanup
  free(recvcounts);
  free(displs);
  free(my_forest_ids);
  free(forest_n_halos);
  free(forest_n_fof);
  free(forest_max_contemp_halos);
  free(forest_max_contemp_fof);
  for (int s = 0; s < n_snaps; s++)
    free(forest_snap_halos[s]);
  free(forest_snap_halos);

  for (int s = 0; s < n_snaps; s++) {
    free(keep[s]);
    free(new_index[s]);
    free(head_target_snap[s]);
    free(head_target_index[s]);
  }
  free(keep);
  free(new_index);
  free(head_target_snap);
  free(head_target_index);
  free(local_n_halos_kept);
  free(local_n_fof_kept);
  free(halo_offset);
  free(global_n_halos);
  free(global_n_fof);

  return EXIT_SUCCESS;
}
