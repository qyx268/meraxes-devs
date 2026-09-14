#include <assert.h>
#include <hdf5_hl.h>
#include <math.h>
#include <string.h>

#include "debug.h"
#include "float_precision_check.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "read_halos.h"
#include "tree_flags.h"
#include "virial_properties.h"

trees_info_t read_trees_info__velociraptor(const int snapshot)
{
  trees_info_t trees_info;

  if (run_globals.mpi_rank == 0) {
    // TODO: This could maybe only ever be done once and stored in run_globals.
    char fname[STRLEN + 34];
    switch (run_globals.params.TreesID) {
      case VELOCIRAPTOR_TREES:
        sprintf(fname, "%s/trees/meraxes_augmented_stats.h5", run_globals.params.SimulationDir);
        break;
      case VELOCIRAPTOR_TREES_AUG:
        sprintf(fname, "%s/augmented_trees/meraxes_augmented_stats.h5", run_globals.params.SimulationDir);
        break;
      default:
        mlog_error("Unrecognised input trees identifier (TreesID).");
        break;
    }

    hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fd < 0) {
      mlog("Failed to open file %s", MLOG_MESG, fname);
      ABORT(EXIT_FAILURE);
    }

    int n_snaps = 0;
    H5LTget_attribute_int(fd, "/", "n_snaps", &n_snaps);
    H5LTget_attribute_int(fd, "/", "n_halos_max", &(trees_info.n_halos_max));
    H5LTget_attribute_int(fd, "/", "n_fof_groups_max", &(trees_info.n_fof_groups_max));

    int* buffer = malloc(sizeof(int) * n_snaps);
    H5LTread_dataset_int(fd, "n_halos", buffer);
    trees_info.n_halos = buffer[snapshot];
    H5LTread_dataset_int(fd, "n_fof_groups", buffer);
    trees_info.n_fof_groups = buffer[snapshot];
    free(buffer);

    H5Fclose(fd);
  }

  // broadcast the snapshot info
  MPI_Bcast(&trees_info, sizeof(trees_info_t), MPI_BYTE, 0, run_globals.mpi_comm);

  return trees_info;
}

static int id_to_ind(long id)
{
  return (int)(((uint64_t)id % (uint64_t)1e12) - 1);
}

static int id_to_snap(long id)
{
  return (int)(id / 1e12l);
}

// N.B. no longer takes/returns a Vvir (neither halo_t nor fof_group_t
// stores it any more, see meraxes.h) or a FOFMvirModifier (mass-ratio-
// modifier support was removed entirely: RequestedMassRatioModifier only
// became 1 if the MassRatioModifier param pointed at a real table file, and
// no simulation config in this repo ever set it, so this was always a
// no-op in practice). Every caller just recomputes Vvir with
// calculate_Vvir(Mvir, Rvir) wherever it's actually needed.
inline static void convert_input_virial_props(float* Mvir,
                                              float* Rvir,
                                              const int len,
                                              const int snapshot,
                                              const bool fof_flag)
{
  double mvir = (double)(*Mvir);
  double rvir = (double)(*Rvir);

  // Update the virial properties for subhalos
  if (mvir == -1) {
    assert(len > 0);
    mvir = calculate_Mvir(mvir, len);
  }

  if (rvir == -1)
    rvir = calculate_Rvir(mvir, snapshot);

  *Mvir = check_float_cast(mvir, fof_flag ? FloatField_FOFGroupMvir : FloatField_HaloMvir);
  *Rvir = check_float_cast(rvir, fof_flag ? FloatField_FOFGroupRvir : FloatField_HaloRvir);
}

// Distributes one column's worth of the chunk rank 0 just read (full_data,
// n_to_read elements of elem_size bytes, meaningful only on rank 0) to their
// owning ranks via MPI_Scatterv, writing the rows this rank owns into
// `dest`. `perm` groups the chunk's rows by destination rank (see the
// caller; meaningful only on rank 0), sendcounts/sdispls are the (element,
// not byte) counts/displacements for that grouping. Every rank, including
// rank 0 itself, receives its own my_recv_count rows into `dest` -- rank 0
// is just an ordinary (non-root-privileged) receiver here as far as the
// receive side goes.
static void scatter_column_from_root(const void* full_data,
                                     const int* perm,
                                     const int* sendcounts,
                                     const int* sdispls,
                                     int mpi_size,
                                     int mpi_rank,
                                     void* send_scratch,
                                     void* dest,
                                     int my_recv_count,
                                     size_t elem_size,
                                     MPI_Comm comm)
{
  int* send_counts_bytes = NULL;
  int* send_displs_bytes = NULL;

  if (mpi_rank == 0) {
    char* send_buf = (char*)send_scratch;
    int total_send = sdispls[mpi_size - 1] + sendcounts[mpi_size - 1];
    for (int k = 0; k < total_send; ++k)
      memcpy(send_buf + (size_t)k * elem_size, (const char*)full_data + (size_t)perm[k] * elem_size, elem_size);

    send_counts_bytes = malloc(sizeof(int) * (size_t)mpi_size);
    send_displs_bytes = malloc(sizeof(int) * (size_t)mpi_size);
    for (int r = 0; r < mpi_size; ++r) {
      send_counts_bytes[r] = sendcounts[r] * (int)elem_size;
      send_displs_bytes[r] = sdispls[r] * (int)elem_size;
    }
  }

  MPI_Scatterv(send_scratch,
               send_counts_bytes,
               send_displs_bytes,
               MPI_BYTE,
               dest,
               my_recv_count * (int)elem_size,
               MPI_BYTE,
               0,
               comm);

  if (mpi_rank == 0) {
    free(send_counts_bytes);
    free(send_displs_bytes);
  }
}

void read_trees__velociraptor(int snapshot,
                              halo_t* halos,
                              int* n_halos,
                              fof_group_t* fof_groups,
                              int* n_fof_groups,
                              int* index_lookup)
{
  switch (run_globals.params.TreesID) {
    case VELOCIRAPTOR_TREES:
      mlog("Reading velociraptor trees for snapshot %d...", MLOG_OPEN, snapshot);
      break;
    case VELOCIRAPTOR_TREES_AUG:
      mlog("Reading velociraptor augmented trees for snapshot %d...", MLOG_OPEN, snapshot);
      break;
    default:
      mlog_error("Unrecognised input trees identifier (TreesID).");
      break;
  }
  
  int n_tree_entries = 0;
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, run_globals.mpi_comm, MPI_INFO_NULL);

  double mass_unit_to_internal = 1.0;
  double scale_factor = -999.;

  *n_halos = 0;
  *n_fof_groups = 0;

  char fname[STRLEN * 2 + 8];
  switch (run_globals.params.TreesID) {
    case VELOCIRAPTOR_TREES:
      sprintf(fname, "%s/trees/%s", run_globals.params.SimulationDir, run_globals.params.CatalogFilePrefix);
      break;
    case VELOCIRAPTOR_TREES_AUG:
      sprintf(fname, "%s/augmented_trees/%s", run_globals.params.SimulationDir, run_globals.params.CatalogFilePrefix);
      break;
    default:
      mlog_error("Unrecognised input trees identifier (TreesID).");
      break;
  }
    
  hid_t fd = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
  if (fd < 0) {
    mlog("Failed to open file %s", MLOG_MESG, fname);
    ABORT(EXIT_FAILURE);
  }
  H5Pclose(plist_id);

  int mpi_size = run_globals.mpi_size;
  int mpi_rank = run_globals.mpi_rank;

  char snap_group_name[9];
  sprintf(snap_group_name, "Snap_%03d", snapshot);
  hid_t snap_group = H5Gopen(fd, snap_group_name, H5P_DEFAULT);

  H5LTget_attribute_int(fd, snap_group_name, "NHalos", &n_tree_entries);

  // check the units
  H5LTget_attribute_double(fd, "Header/Units", "Mass_unit_to_solarmass", &mass_unit_to_internal);
  mass_unit_to_internal /= 1.0e10;
  H5LTget_attribute_double(fd, snap_group_name, "scalefactor", &scale_factor);

  // Currently the chunk size is 10000, improving ~10% w.r.t. no chunk or 1k chunk
  int buffer_size = (n_tree_entries > 100000) ? n_tree_entries / 10 : 10000;
  buffer_size = buffer_size > n_tree_entries ? n_tree_entries : buffer_size;

  long* Head = malloc(sizeof(long) * buffer_size);
  long* hostHaloID = malloc(sizeof(long) * buffer_size);
  float* Mass_200crit = malloc(sizeof(float) * buffer_size);
  float* Mass_tot = malloc(sizeof(float) * buffer_size);
  float* R_200crit = malloc(sizeof(float) * buffer_size);
  float* Vmax = malloc(sizeof(float) * buffer_size);
  float* Xc = malloc(sizeof(float) * buffer_size);
  float* Yc = malloc(sizeof(float) * buffer_size);
  float* Zc = malloc(sizeof(float) * buffer_size);
  // N.B. VXc/VYc/VZc are no longer read -- neither halo_t nor galaxy_t
  // tracks a Vel field any more (it was only ever used by the output
  // writer, which no longer has anywhere to source it from).
  float* AngMom = malloc(sizeof(float) * buffer_size);
  unsigned long* ID = malloc(sizeof(unsigned long) * buffer_size); 
  unsigned long* npart = malloc(sizeof(unsigned long) * buffer_size);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); // or H5FD_MPIO_COLLECTIVE?
  hid_t fspace_id = H5Screate_simple(1, (hsize_t[1]){ n_tree_entries }, NULL);

  double hubble_h = run_globals.params.Hubble_h;
  double box_size = run_globals.params.BoxSize;

  // Rank 0 always does the actual disk I/O (one big sequential read per
  // column per chunk, same access pattern as before any of this
  // parallelization work) and then targets its distribution using the
  // global forest-id -> owning-rank map, instead of broadcasting the full
  // chunk to everyone. All of this bookkeeping is only ever touched on rank
  // 0 -- every other rank just receives its own share via MPI_Scatterv.
  bool have_forest_map = (run_globals.RequestedForestId != NULL);
  int* sendcounts = NULL;
  int* sdispls = NULL;
  int* dest_rank = NULL;
  int* perm = NULL;
  void* send_scratch = NULL;
  if (mpi_rank == 0) {
    sendcounts = malloc(sizeof(int) * mpi_size);
    sdispls = malloc(sizeof(int) * mpi_size);
    dest_rank = malloc(sizeof(int) * buffer_size);
    perm = malloc(sizeof(int) * buffer_size);
    // Big enough for any column (long/unsigned long are the largest of the
    // three types we read).
    send_scratch = malloc(sizeof(long) * buffer_size);
  }

  int n_read = 0;
  int n_to_read = buffer_size;
  while (n_read < n_tree_entries) {
    int n_remaining = n_tree_entries - n_read;
    if (n_remaining < n_to_read)
      n_to_read = n_remaining;


    // select a hyperslab in the filespace
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, (hsize_t[1]){ n_read }, NULL, (hsize_t[1]){ n_to_read }, NULL);
    hid_t memspace_id = H5Screate_simple(1, (hsize_t[1]){ n_to_read }, NULL);

    // Rank 0 reads ForestID first -- everything else (who owns each row, so
    // how to group this chunk for MPI_Scatterv) depends on it -- and
    // computes the per-row destination grouping. Nothing here is needed by
    // any other rank, so it's all local to rank 0.
    if (mpi_rank == 0) {
      long* full_forestid = malloc(sizeof(long) * (size_t)n_to_read);
      {
        hid_t dset_id = H5Dopen(snap_group, "ForestID", H5P_DEFAULT);
        herr_t status = H5Dread(dset_id, H5T_NATIVE_LONG, memspace_id, fspace_id, plist_id, full_forestid);
        assert(status >= 0);
        H5Dclose(dset_id);
      }

      for (int r = 0; r < mpi_size; ++r)
        sendcounts[r] = 0;
      for (int ii = 0; ii < n_to_read; ++ii) {
        dest_rank[ii] = have_forest_map ? get_forest_owner_rank(full_forestid[ii]) : 0;
        if (dest_rank[ii] >= 0)
          sendcounts[dest_rank[ii]]++;
      }
      free(full_forestid);

      {
        int offset = 0;
        for (int r = 0; r < mpi_size; ++r) {
          sdispls[r] = offset;
          offset += sendcounts[r];
        }
      }
      {
        // fill_pos tracks, per destination rank, the next free slot in perm
        // -- starts at sdispls and advances as each row of that destination
        // is placed, so rows destined for the same rank stay in their
        // original (original file) order.
        int* fill_pos = malloc(sizeof(int) * mpi_size);
        memcpy(fill_pos, sdispls, sizeof(int) * mpi_size);
        for (int ii = 0; ii < n_to_read; ++ii)
          if (dest_rank[ii] >= 0)
            perm[fill_pos[dest_rank[ii]]++] = ii;
        free(fill_pos);
      }
    }

    // Every rank (including rank 0) learns how many rows it's about to
    // receive this chunk before the real Scatterv calls, since each rank
    // has to supply its own receive count up front.
    int my_recv_count = 0;
    MPI_Scatter(sendcounts, 1, MPI_INT, &my_recv_count, 1, MPI_INT, 0, run_globals.mpi_comm);

#define READ_AND_SCATTER(name, type, h5type)                                                                          \
{                                                                                                                      \
  type* full_col = NULL;                                                                                              \
  if (mpi_rank == 0) {                                                                                                 \
    full_col = malloc(sizeof(type) * (size_t)n_to_read);                                                              \
    hid_t dset_id = H5Dopen(snap_group, #name, H5P_DEFAULT);                                                           \
    herr_t status = H5Dread(dset_id, h5type, memspace_id, fspace_id, plist_id, full_col);                              \
    assert(status >= 0);                                                                                               \
    H5Dclose(dset_id);                                                                                                 \
  }                                                                                                                    \
  scatter_column_from_root(full_col,                                                                                   \
                           perm,                                                                                       \
                           sendcounts,                                                                                 \
                           sdispls,                                                                                    \
                           mpi_size,                                                                                   \
                           mpi_rank,                                                                                   \
                           send_scratch,                                                                               \
                           name,                                                                                       \
                           my_recv_count,                                                                              \
                           sizeof(type),                                                                               \
                           run_globals.mpi_comm);                                                                      \
  if (mpi_rank == 0)                                                                                                   \
    free(full_col);                                                                                                    \
}                                                                                                                      \

    READ_AND_SCATTER(Head,         long, H5T_NATIVE_LONG);
    READ_AND_SCATTER(hostHaloID,   long, H5T_NATIVE_LONG);
    READ_AND_SCATTER(Mass_200crit, float, H5T_NATIVE_FLOAT);
    READ_AND_SCATTER(Mass_tot,     float, H5T_NATIVE_FLOAT);
    READ_AND_SCATTER(R_200crit,    float, H5T_NATIVE_FLOAT);
    READ_AND_SCATTER(Vmax,         float, H5T_NATIVE_FLOAT);
    READ_AND_SCATTER(Xc,           float, H5T_NATIVE_FLOAT);
    READ_AND_SCATTER(Yc,           float, H5T_NATIVE_FLOAT);
    READ_AND_SCATTER(Zc,           float, H5T_NATIVE_FLOAT);
    READ_AND_SCATTER(AngMom,       float, H5T_NATIVE_FLOAT);
    READ_AND_SCATTER(ID,           unsigned long, H5T_NATIVE_ULONG);
    READ_AND_SCATTER(npart,        unsigned long, H5T_NATIVE_ULONG);

    H5Sclose(memspace_id);

    // Every row here was already distributed to this rank via
    // MPI_Scatterv above because it owns that row's forest (see
    // get_forest_owner_rank()), so there's no per-row forest-membership
    // filter needed any more -- every row in [0, my_recv_count) is kept.
    for (int ii = 0; ii < my_recv_count; ++ii) {
      {
        halo_t* halo = &(halos[*n_halos]);

        // N.B. halo_t no longer stores the raw catalogue ID (it was only
        // ever used to seed gal->ID/galout->HaloID, both now dropped). The
        // local `ID` array read above is still needed below though, for the
        // self-referencing "Head[ii]==ID[ii]" terminal-halo check.
        halo->DescIndex = id_to_ind(Head[ii]);

        // N.B. halo_t no longer stores ProgIndex -- see meraxes.h for why
        // (write-only: FlagIgnoreProgIndex is 1 for every simulation config
        // in this repo, and the alternative computation below was never
        // implemented anyway, since `Tail` is never read from the file).
        //if (!run_globals.params.FlagIgnoreProgIndex)
        //  halo->ProgIndex = id_to_ind(Tail[ii]);

        halo->NextHaloInFOFGroup = NULL;

        // Any other tree flags need to be set using both the current and
        // progenitor halo information (stored in the galaxy), therefore we
        // need to leave setting those until later...
        if (run_globals.params.FlagIgnoreProgIndex)
          halo->TreeFlags = TREE_CASE_NO_PROGENITORS;
        //else
        //  halo->TreeFlags = (unsigned long)Tail[ii] != ID[ii] ? 0 : TREE_CASE_NO_PROGENITORS;

        // N.B. Type/SnapOffset must be set after the TreeFlags assignment(s)
        // above, since they're packed into TreeFlags's own high bits and a
        // plain `=` (as opposed to `|=`) would otherwise wipe them out.
        halo_set_type(halo, hostHaloID[ii] == -1 ? 0 : 1);
        halo_set_snap_offset(halo, id_to_snap(Head[ii]) - snapshot);

        // Here we have a cyclic pointer, indicating that this halo's life ends here
        if ((unsigned long)Head[ii] == ID[ii])
          halo->DescIndex = -1;

        // N.B. no longer "ii + n_read" -- after distribution, ii no longer
        // indexes contiguously into this chunk's original file rows (rows
        // were regrouped by owning rank), so the original row index has to
        // come from this halo's own ID instead (same id_to_ind() convention
        // used above for Head/hostHaloID; ID is distributed alongside
        // everything else, so it's still correct post-distribution).
        if (index_lookup)
          index_lookup[*n_halos] = id_to_ind((long)ID[ii]);

        // TODO: What masses and radii should I use for centrals (inclusive vs. exclusive etc.)?
        if (halo_get_type(halo) == 0) {
          fof_group_t* fof_group = &fof_groups[*n_fof_groups];
          
          // This check is to ensure sensible values of mass_200crit and avoid having 
          // very weird halos. Put this check back if you feel that the N-body is weird.
          
          /*if ((Mass_200crit[ii] < 5 * Mass_tot[ii])) {
              fof_group->Mvir = Mass_200crit[ii] * hubble_h * mass_unit_to_internal;;
              fof_group->Rvir = R_200crit[ii] * hubble_h;
          }*/
          //else {
            // BELOW_VIRIAL_THRESHOLD merger halo swammping
          if (Mass_200crit[ii] <= 0) {
            halo->TreeFlags |= TREE_CASE_BELOW_VIRIAL_THRESHOLD;
            fof_group->Mvir = check_float_cast((double)Mass_tot[ii] * hubble_h * mass_unit_to_internal, FloatField_FOFGroupMvir);
            fof_group->Rvir = check_float_cast(-1, FloatField_FOFGroupRvir);
          }

          else {
            fof_group->Mvir = check_float_cast((double)Mass_200crit[ii] * hubble_h * mass_unit_to_internal, FloatField_FOFGroupMvir);
            fof_group->Rvir = check_float_cast((double)R_200crit[ii] * hubble_h, FloatField_FOFGroupRvir);
          }

          convert_input_virial_props(&fof_group->Mvir, &fof_group->Rvir, -1, snapshot, true);

          halo->FOFGroup = &(fof_groups[*n_fof_groups]);
          fof_groups[(*n_fof_groups)++].FirstHalo = halo;
        } else {
          // We can take advantage of the fact that host halos always
          // seem to appear before their subhalos (checked below) in the
          // trees to immediately connect FOF group members.
          int host_index = id_to_ind(hostHaloID[ii]);

          if (index_lookup)
            host_index = find_original_index(host_index, index_lookup, *n_halos);

          assert(host_index > -1);
          assert(host_index < *n_halos);

          halo_t* prev_halo = &halos[host_index];
          halo->FOFGroup = prev_halo->FOFGroup;

          while (prev_halo->NextHaloInFOFGroup != NULL)
            prev_halo = prev_halo->NextHaloInFOFGroup;

          prev_halo->NextHaloInFOFGroup = halo;
        }

        halo->Len = (int)npart[ii];
        halo->Pos[0] = apply_pbc_pos(Xc[ii] * hubble_h / scale_factor);
        halo->Pos[1] = apply_pbc_pos(Yc[ii] * hubble_h / scale_factor);
        halo->Pos[2] = apply_pbc_pos(Zc[ii] * hubble_h / scale_factor);
        halo->Vmax = Vmax[ii];

        // TODO: What masses and radii should I use for satellites (inclusive vs. exclusive etc.)?
        halo->Mvir = check_float_cast((double)Mass_tot[ii] * hubble_h * mass_unit_to_internal, FloatField_HaloMvir);
        halo->Rvir = check_float_cast(-1, FloatField_HaloRvir);
        convert_input_virial_props(&halo->Mvir, &halo->Rvir, -1, snapshot, false);

        halo->AngMom = AngMom[ii] * hubble_h;
        halo->Galaxy = NULL;

        (*n_halos)++;
      }
    }

    n_read += n_to_read;
  }

  free(Head);
  free(hostHaloID);
  free(Mass_200crit);
  free(Mass_tot);
  free(R_200crit);
  free(Vmax);
  free(Xc);
  free(Yc);
  free(Zc);
  free(AngMom);
  free(ID);
  free(npart);
  if (mpi_rank == 0) {
    free(sendcounts);
    free(sdispls);
    free(dest_rank);
    free(perm);
    free(send_scratch);
  }
  H5Pclose(plist_id);
  H5Sclose(fspace_id);
  H5Gclose(snap_group);
  H5Fclose(fd);

  mlog("...done", MLOG_CLOSE);
}
