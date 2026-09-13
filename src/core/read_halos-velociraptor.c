#include <assert.h>
#include <hdf5_hl.h>
#include <math.h>

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

  long* ForestID = malloc(sizeof(long) * buffer_size);
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
  // Back to independent, not collective: every rank still reads its own
  // genuinely disjoint slice of each column concurrently (see below), which
  // is the actual parallelism win. But H5FD_MPIO_COLLECTIVE's two-phase I/O
  // aggregation is only a net win when tuned to the filesystem (Lustre
  // striping/aggregator-count hints), and H5Pset_fapl_mpio above passes no
  // such hints (MPI_INFO_NULL) -- untuned, that machinery is pure
  // coordination overhead on top of the read, and was observed to make
  // things slower rather than faster. Independent transfers don't pay that
  // coordination cost; each rank's H5Dread just proceeds on its own.
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  hid_t fspace_id = H5Screate_simple(1, (hsize_t[1]){ n_tree_entries }, NULL);

  double hubble_h = run_globals.params.Hubble_h;
  double box_size = run_globals.params.BoxSize;

  // Even split of each chunk's rows across ranks, recomputed once per chunk
  // and reused for every column's hyperslab selection and MPI_Allgatherv.
  int* rank_n_rows = malloc(sizeof(int) * mpi_size);
  int* rank_row_offset = malloc(sizeof(int) * mpi_size);
  int* recvcounts_bytes = malloc(sizeof(int) * mpi_size);
  int* displs_bytes = malloc(sizeof(int) * mpi_size);

  int n_read = 0;
  int n_to_read = buffer_size;
  while (n_read < n_tree_entries) {
    int n_remaining = n_tree_entries - n_read;
    if (n_remaining < n_to_read)
      n_to_read = n_remaining;

    {
      int base = n_to_read / mpi_size;
      int rem = n_to_read % mpi_size;
      int offset = 0;
      for (int r = 0; r < mpi_size; ++r) {
        rank_n_rows[r] = base + (r < rem ? 1 : 0);
        rank_row_offset[r] = offset;
        offset += rank_n_rows[r];
      }
    }
    int my_n_rows = rank_n_rows[mpi_rank];
    int my_row_offset = rank_row_offset[mpi_rank];

    // select this rank's own (disjoint) hyperslab of the current chunk
    H5Sselect_hyperslab(
      fspace_id, H5S_SELECT_SET, (hsize_t[1]){ n_read + my_row_offset }, NULL, (hsize_t[1]){ my_n_rows }, NULL);
    hid_t memspace_id = H5Screate_simple(1, (hsize_t[1]){ my_n_rows }, NULL);

#define READ_TREE_ENTRY_PROP(name, type, h5type)                                                                       \
{                                                                                                                      \
  type* my_slice = malloc(sizeof(type) * (size_t)my_n_rows);                                                          \
  hid_t dset_id = H5Dopen(snap_group, #name, H5P_DEFAULT);                                                            \
  herr_t status = H5Dread(dset_id, h5type, memspace_id, fspace_id, plist_id, my_slice);                                \
  assert(status >= 0);                                                                                                 \
  H5Dclose(dset_id);                                                                                                   \
  for (int r = 0; r < mpi_size; ++r) {                                                                                 \
    recvcounts_bytes[r] = rank_n_rows[r] * (int)sizeof(type);                                                          \
    displs_bytes[r] = rank_row_offset[r] * (int)sizeof(type);                                                          \
  }                                                                                                                     \
  MPI_Allgatherv(my_slice,                                                                                             \
                 my_n_rows * (int)sizeof(type),                                                                        \
                 MPI_BYTE,                                                                                             \
                 name,                                                                                                 \
                 recvcounts_bytes,                                                                                     \
                 displs_bytes,                                                                                         \
                 MPI_BYTE,                                                                                             \
                 run_globals.mpi_comm);                                                                                \
  free(my_slice);                                                                                                      \
}                                                                                                                      \

    READ_TREE_ENTRY_PROP(ForestID,     long, H5T_NATIVE_LONG);
    READ_TREE_ENTRY_PROP(Head,         long, H5T_NATIVE_LONG);
    READ_TREE_ENTRY_PROP(hostHaloID,   long, H5T_NATIVE_LONG);
    READ_TREE_ENTRY_PROP(Mass_200crit, float, H5T_NATIVE_FLOAT);
    READ_TREE_ENTRY_PROP(Mass_tot,     float, H5T_NATIVE_FLOAT);
    READ_TREE_ENTRY_PROP(R_200crit,    float, H5T_NATIVE_FLOAT);
    READ_TREE_ENTRY_PROP(Vmax,         float, H5T_NATIVE_FLOAT);
    READ_TREE_ENTRY_PROP(Xc,           float, H5T_NATIVE_FLOAT);
    READ_TREE_ENTRY_PROP(Yc,           float, H5T_NATIVE_FLOAT);
    READ_TREE_ENTRY_PROP(Zc,           float, H5T_NATIVE_FLOAT);
    READ_TREE_ENTRY_PROP(AngMom,       float, H5T_NATIVE_FLOAT);
    READ_TREE_ENTRY_PROP(ID,           unsigned long, H5T_NATIVE_ULONG);
    READ_TREE_ENTRY_PROP(npart,        unsigned long, H5T_NATIVE_ULONG);

    H5Sclose(memspace_id);

    for (int ii = 0; ii < n_to_read; ++ii) {
      bool keep_this_halo = true;

      if ((run_globals.RequestedForestId != NULL) && (bsearch(&(ForestID[ii]),
                                                              run_globals.RequestedForestId,
                                                              (size_t)run_globals.NRequestedForests,
                                                              sizeof(long),
                                                              compare_longs)) == NULL)
        keep_this_halo = false;

      if (keep_this_halo) {
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

        if (index_lookup)
          index_lookup[*n_halos] = ii + n_read;

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

  free(ForestID);
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
  free(rank_n_rows);
  free(rank_row_offset);
  free(recvcounts_bytes);
  free(displs_bytes);
  H5Pclose(plist_id);
  H5Sclose(fspace_id);
  H5Gclose(snap_group);
  H5Fclose(fd);

  mlog("...done", MLOG_CLOSE);
}
