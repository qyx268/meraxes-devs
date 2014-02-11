#include "meraxes.h"
#include "tree_flags.h"
#include <hdf5.h>
#include <hdf5_hl.h>

int get_corrected_snapshot(run_globals_t *run_globals, int snapshot)
{
  int  total_sim_snaps = run_globals->params.TotalSimSnaps;
  int  n_every_snaps   = run_globals->params.NEverySnap;

  // Calculate the actual (unsampled) simulation snapshot.
  if (n_every_snaps>1)
    return total_sim_snaps - ((int)((total_sim_snaps+1)/n_every_snaps))*n_every_snaps + snapshot*n_every_snaps;
  else
    return snapshot;
}

static inline bool is_ghost(int flags)
{
 if ((flags & TREE_CASE_GHOST)==TREE_CASE_GHOST)
   return true;
 else
   return false;
}

static void halo_catalog_filename(
  char *simulation_dir,
  char *catalog_file_prefix,
  int   snapshot,
  char *group_type,
  int   sub,
  int  *i_layout,
  char *fname)
{

  bool flag_success = false;
  FILE *fin;

  // if we need to determine the filename structure...
  if (*i_layout==-1)
  {
    for (*i_layout=0; (*i_layout<4) && (flag_success==false); (*i_layout)++)
    {
      if (*i_layout==0)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type, sub);
      else if (*i_layout==1)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type);
      else if (*i_layout==2)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, sub);
      else if (*i_layout==3)
        sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type);

      if ((fin = fopen(fname, "rb"))!=NULL)
      {
        flag_success = true;
        fclose(fin);
        break;
      }
    }
  }

  // ensure we have a valid i_layout value.
  if (*i_layout<0 && *i_layout>3)
  {
    fprintf(stderr, "cannot resolve catalogue filename.\n");
    ABORT(EXIT_FAILURE);
  }

  // provide the correct filename
  if (*i_layout==0)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type, sub);
  else if (*i_layout==1)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type, catalog_file_prefix, snapshot, group_type);
  else if (*i_layout==2)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties.%d", simulation_dir, catalog_file_prefix, snapshot, group_type, sub);
  else if (*i_layout==3)
    sprintf(fname, "%s/catalogs/%s_%03d.catalog_%s_properties", simulation_dir, catalog_file_prefix, snapshot, group_type);

}


static void inline read_catalogs_header(
  FILE *fin,
  int  *i_file,
  int  *N_files,
  int  *N_halos_file,
  int  *N_halos_total )
{
  fread(i_file       , sizeof(int), 1, fin);
  fread(N_files      , sizeof(int), 1, fin);
  fread(N_halos_file , sizeof(int), 1, fin);
  fread(N_halos_total, sizeof(int), 1, fin);
}


static void read_catalog_halos(
    FILE           **fin,
    char           *simulation_dir,
    char           *catalog_file_prefix,
    int             snapshot,
    int            *flayout_switch,
    int            *i_file,
    int            *N_halos_file,
    int            *i_halo,
    catalog_halo_t *halo,
    int             N_to_read)
{

  char                 fname[STRLEN];
  int                  dummy;
  int                  N_from_this_file;

  // Is this the first read?
  if((*fin)==NULL)
  {
    halo_catalog_filename(simulation_dir, catalog_file_prefix, snapshot, "subgroups", *i_file, flayout_switch, fname);
    *fin = fopen(fname, "rb");
    if (*fin==NULL)
    {
      SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
      ABORT(34494);
    }
    read_catalogs_header(*fin, &dummy, &dummy, N_halos_file, &dummy);
  }

  // Have we already read all the halos in this file?
  if((*i_halo)>=(*N_halos_file))
  {
    // SID_log("***i_halo = %d, N_halos_file = %d", SID_LOG_COMMENT, (*i_halo), (*N_halos_file));
    fclose(*fin);
    (*i_file)++;
    (*i_halo) = 0;
    halo_catalog_filename(simulation_dir, catalog_file_prefix, snapshot, "subgroups", *i_file, flayout_switch, fname);
    *fin = fopen(fname, "rb");
    if (*fin==NULL)
    {
      SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
      ABORT(34494);
    }
    read_catalogs_header(*fin, &dummy, &dummy, N_halos_file, &dummy);
  }

  // Read in as many halos as we can from this file
  if((*i_halo + N_to_read) <= *N_halos_file)
  {
    fread(halo, sizeof(catalog_halo_t), N_to_read, *fin);
    *i_halo += N_to_read;
  }
  else
  {
    // read in as many as we can from this file and then get the rest from the next file
    N_from_this_file = (*N_halos_file)- *i_halo;

    // DEBUG
    // SID_log("Spilling over to next file: i_halo = %d, N_to_read = %d,"\
    //     "N_halos_file = %d, N_from_this_file = %d", SID_LOG_COMMENT, *i_halo,
    //     N_to_read, *N_halos_file, N_from_this_file);

    fread(halo, sizeof(catalog_halo_t), N_from_this_file, *fin);
    *i_halo += N_from_this_file;
    N_to_read -= N_from_this_file;
    read_catalog_halos(fin, simulation_dir, catalog_file_prefix, snapshot, flayout_switch, i_file, N_halos_file, i_halo, halo, N_to_read);
  }

}


static void inline convert_input_halo_units(run_globals_t *run_globals, halo_t *halo, int snapshot)
{
  halo->Mvir /= 1.0e10;

  // Update the virial properties
  halo->Mvir = calculate_Mvir(run_globals, halo);
  halo->Rvir = calculate_Rvir(run_globals, halo, halo->Mvir, snapshot);
  halo->Vvir = calculate_Vvir(run_globals, halo->Mvir, halo->Rvir);
}



int compare_ints(const void *a, const void *b)
{
  return *((int *)a) - *((int *)b);
}

//! Buffered read of hdf5 trees into halo structures
static void read_trees_and_catalogs(
  run_globals_t *run_globals,
  int            snapshot,
  hid_t          fd,
  halo_t        *halo,
  int            N_halos,
  fof_group_t   *fof_group,
  int            N_fof_groups,
  int           *requested_forest_id,
  int            N_requested_forests,
  int           *N_halos_kept,
  int           *N_fof_groups_kept)
{
  // I guess this should ideally be equal to the chunk size of the input hdf5 file...
  int buffer_size = 1000;
  int N_read = 0;
  int N_to_read = 0;
  bool keep_flag;

  FILE *fin_catalogs = NULL;
  int flayout_switch = -1;
  int i_catalog_file = 0;
  int N_halos_in_catalog_file = 0;
  int i_halo_in_catalog_file = 0;
  char catalog_file_prefix[50];
  char simulation_dir[STRLEN];
  catalog_halo_t *catalog_buffer;

  sprintf(catalog_file_prefix, "%s", run_globals->params.CatalogFilePrefix);
  sprintf(simulation_dir, "%s", run_globals->params.SimulationDir);

  tree_entry_t *tree_buffer;
  tree_buffer = SID_malloc(sizeof(tree_entry_t) * buffer_size);
  catalog_buffer = SID_malloc(sizeof(catalog_halo_t) * buffer_size);

  *N_halos_kept = 0;
  *N_fof_groups_kept = 0;

  // DEBUG
  // SID_log("Calling read_trees_and_catalogs() with:", SID_LOG_OPEN);
  // SID_log("snapshot = %d", SID_LOG_COMMENT, snapshot);
  // SID_log("N_halos = %d", SID_LOG_COMMENT, N_halos);
  // SID_log("N_fof_groups = %d", SID_LOG_COMMENT, N_fof_groups);
  // SID_log("N_requested_forests = %d", SID_LOG_COMMENT, N_requested_forests);
  // SID_log("---", SID_LOG_CLOSE);

  size_t dst_size = sizeof(tree_entry_t);
  size_t dst_offsets[9] = {
    HOFFSET(tree_entry_t, id),
    HOFFSET(tree_entry_t, flags),
    HOFFSET(tree_entry_t, desc_id),
    HOFFSET(tree_entry_t, tree_id),
    HOFFSET(tree_entry_t, file_offset),
    HOFFSET(tree_entry_t, desc_index),
    HOFFSET(tree_entry_t, central_index),
    HOFFSET(tree_entry_t, forest_id),
    HOFFSET(tree_entry_t, fof_mvir) };
  size_t dst_sizes[9] = {
    sizeof(tree_buffer[0].id),
    sizeof(tree_buffer[0].flags),
    sizeof(tree_buffer[0].desc_id),
    sizeof(tree_buffer[0].tree_id),
    sizeof(tree_buffer[0].file_offset),
    sizeof(tree_buffer[0].desc_index),
    sizeof(tree_buffer[0].central_index),
    sizeof(tree_buffer[0].forest_id),
    sizeof(tree_buffer[0].fof_mvir) };


  keep_flag = true;
  while(N_read<N_halos)
  {

    if((N_halos-N_read) >= buffer_size)
      N_to_read = buffer_size;
    else
      N_to_read = N_halos-N_read;

    // read in a tree_buffer of the trees
    H5TBread_records(fd, "trees", N_read, (hsize_t)N_to_read, dst_size, dst_offsets, dst_sizes, tree_buffer);

    // read in the corresponding catalog entrys
    read_catalog_halos(&fin_catalogs, simulation_dir, catalog_file_prefix,
        snapshot, &flayout_switch, &i_catalog_file, &N_halos_in_catalog_file,
        &i_halo_in_catalog_file, catalog_buffer, N_to_read);

    // paste the data into the halo structures
    for(int jj=0; jj<N_to_read; jj++)
    {

      if(requested_forest_id!=NULL)
      {
        if(bsearch(&(tree_buffer[jj].forest_id), requested_forest_id,
              (size_t)N_requested_forests, sizeof(int), compare_ints) != NULL)
          keep_flag = true;
        else
          keep_flag = false;
      }

      if(keep_flag)
      {
        halo[*N_halos_kept].ID = tree_buffer[jj].id;
        halo[*N_halos_kept].TreeFlags = tree_buffer[jj].flags;
        halo[*N_halos_kept].SnapOffset = tree_buffer[jj].file_offset;
        halo[*N_halos_kept].DescIndex = tree_buffer[jj].desc_index;
        halo[*N_halos_kept].Mvir = tree_buffer[jj].fof_mvir;  // this will be overwritten for type>0 halos later
        halo[*N_halos_kept].NextHaloInFOFGroup = NULL;

        if((N_read+jj) == tree_buffer[jj].central_index)
        {
          halo[*N_halos_kept].Type = 0;
          fof_group[(*N_fof_groups_kept)++].FirstHalo = &(halo[*N_halos_kept]);
        }
        else
        {
          halo[*N_halos_kept].Type = 1;
          halo[*N_halos_kept-1].NextHaloInFOFGroup = &(halo[*N_halos_kept]);
        }

        halo[*N_halos_kept].FOFGroup = &(fof_group[(*N_fof_groups_kept)-1]);

        // paste in the halo properties
        halo[*N_halos_kept].id_MBP             = catalog_buffer[jj].id_MBP;
        halo[*N_halos_kept].Len                = catalog_buffer[jj].n_particles;
        halo[*N_halos_kept].Pos[0]             = catalog_buffer[jj].position_MBP[0];
        halo[*N_halos_kept].Pos[1]             = catalog_buffer[jj].position_MBP[1];
        halo[*N_halos_kept].Pos[2]             = catalog_buffer[jj].position_MBP[2];
        halo[*N_halos_kept].Vel[0]             = catalog_buffer[jj].velocity_COM[0];
        halo[*N_halos_kept].Vel[1]             = catalog_buffer[jj].velocity_COM[1];
        halo[*N_halos_kept].Vel[2]             = catalog_buffer[jj].velocity_COM[2];
        halo[*N_halos_kept].Rvir               = catalog_buffer[jj].R_vir;
        halo[*N_halos_kept].Rhalo              = catalog_buffer[jj].R_halo;
        halo[*N_halos_kept].Rmax               = catalog_buffer[jj].R_max;
        halo[*N_halos_kept].Vmax               = catalog_buffer[jj].V_max;
        halo[*N_halos_kept].VelDisp            = catalog_buffer[jj].sigma_v;
        halo[*N_halos_kept].Spin[0]            = catalog_buffer[jj].spin[0];
        halo[*N_halos_kept].Spin[1]            = catalog_buffer[jj].spin[1];
        halo[*N_halos_kept].Spin[2]            = catalog_buffer[jj].spin[2];
        halo[*N_halos_kept].Galaxy             = NULL;
        if(halo[*N_halos_kept].Type > 0)
        {
          halo[*N_halos_kept].Mvir             = catalog_buffer[jj].M_vir;
        }

        convert_input_halo_units(run_globals, &(halo[*N_halos_kept]), snapshot);

        (*N_halos_kept)++;
      }

    }

    N_read += N_to_read;
  }

  // close the catalogs file
  if(fin_catalogs)
    fclose(fin_catalogs);

  // free the buffers
  SID_free(SID_FARG catalog_buffer);
  SID_free(SID_FARG tree_buffer);

}


static trees_info_t read_trees_info(hid_t fd)
{

  trees_info_t trees_info;

  H5LTget_attribute_int(fd, "trees", "n_step", &(trees_info.n_step));
  H5LTget_attribute_int(fd, "trees", "n_search", &(trees_info.n_search));
  H5LTget_attribute_int(fd, "trees", "n_halos", &(trees_info.n_halos));
  H5LTget_attribute_int(fd, "trees", "n_halos_max", &(trees_info.n_halos_max));
  H5LTget_attribute_int(fd, "trees", "max_tree_id", &(trees_info.max_tree_id));
  H5LTget_attribute_int(fd, "trees", "n_fof_groups", &(trees_info.n_fof_groups));
  H5LTget_attribute_int(fd, "trees", "n_fof_groups_max", &(trees_info.n_fof_groups_max));

  return trees_info;

}


static void read_forests_info(run_globals_t *run_globals, int *requested_forest_id, int N_requested_forests)
{

  char fname[STRLEN];
  hid_t fin;
  int n_forests;
  int *forest_id;
  int *max_contemp_halo;
  int *max_contemp_fof;
  int max_halos;
  int max_fof_groups;

  sprintf(fname, "%s/trees/forests_info.hdf5", run_globals->params.SimulationDir);
  if ((fin = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(EXIT_FAILURE);
  }

  // find out how many forests there are
  H5LTget_attribute_int(fin, "info", "n_forests", &n_forests);

  // allocate the arrays
  forest_id = SID_malloc(sizeof(int) * n_forests);
  max_contemp_halo = SID_malloc(sizeof(int) * n_forests);
  max_contemp_fof = SID_malloc(sizeof(int) * n_forests);

  // read in the max number of contemporaneous halos and groups and the forest ids
  H5LTread_dataset_int(fin, "info/max_contemporaneous_halos", max_contemp_halo);
  H5LTread_dataset_int(fin, "info/max_contemporaneous_fof_groups", max_contemp_fof);
  H5LTread_dataset_int(fin, "info/forest_id", forest_id);

  // close the file
  H5Fclose(fin);

  // loop through and tot up the max number of halos and fof_groups we will need to allocate
  max_halos = 0;
  max_fof_groups = 0;
  for(int i_forest=0, i_req=0; (i_forest<n_forests) && (i_req<N_requested_forests); i_forest++)
  {
    if(forest_id[i_forest] == requested_forest_id[i_req])
    {
      max_halos += max_contemp_halo[i_forest];
      max_fof_groups += max_contemp_fof[i_forest];
      i_req++;
    }
  }

  // store the maximum number of halos and fof groups needed at any one snapshot
  run_globals->N_halos_max = max_halos;
  run_globals->N_fof_groups_max = max_fof_groups;

  // free the arrays
  SID_free(SID_FARG forest_id);
  SID_free(SID_FARG max_contemp_fof);
  SID_free(SID_FARG max_contemp_halo);

}


trees_info_t read_halos(
  run_globals_t  *run_globals,
  int                  snapshot,
  halo_t        **halo,
  fof_group_t   **fof_group)
{

  int             N_halos;                 //!< Number of halos
  char            fname[STRLEN];
  int             corrected_snapshot;
  int             N_fof_groups;
  trees_info_t    trees_info;
  hid_t           fin_trees;

  int             N_halos_kept;
  int             N_fof_groups_kept;
  bool            subsample_trees = true;  // TEMPORARY FOR DEVELOPMENT
  int             N_requested_forests = 1;  // TEMPORARY FOR DEVELOPMENT
  int             requested_forest_id[1] = {243664};  // TEMPORARY FOR DEVELOPMENT

  SID_log("Reading snapshot %d (z=%.2f) trees and halos...", SID_LOG_OPEN|SID_LOG_TIMER, snapshot, run_globals->ZZ[snapshot]);

  corrected_snapshot = get_corrected_snapshot(run_globals, snapshot);

  // open the tree file
  sprintf(fname, "%s/trees/horizontal_trees_%03d.hdf5", run_globals->params.SimulationDir, corrected_snapshot);
  if ((fin_trees = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
  {
    SID_log("Failed to open file %s", SID_LOG_COMMENT, fname);
    ABORT(EXIT_FAILURE);
  }

  // read the info attributes
  trees_info = read_trees_info(fin_trees);
  N_halos = trees_info.n_halos;
  N_fof_groups = trees_info.n_fof_groups;

  // If necessary, allocate the halo array
  if(*halo == NULL)
  {
    // if required, read the forest info and calculate the maximum number of halos and fof groups
    if(subsample_trees)
      read_forests_info(run_globals, requested_forest_id, N_requested_forests);
    else
    {
      run_globals->N_halos_max = trees_info.n_halos_max;
      run_globals->N_fof_groups_max = trees_info.n_fof_groups_max;
    }

    SID_log("Allocating halo array with %d elements...", SID_LOG_COMMENT, run_globals->N_halos_max);
    *halo = SID_malloc(sizeof(halo_t) * run_globals->N_halos_max);
  }

  // Allocate the fof_group array if necessary
  if(*fof_group == NULL)
  {
    SID_log("Allocating fof_group array with %d elements...", SID_LOG_COMMENT, run_globals->N_fof_groups_max);
    *fof_group = SID_malloc(sizeof(fof_group_t) * run_globals->N_fof_groups_max);
  }

  // reset the fof group pointers
  for(int ii=0; ii<run_globals->N_fof_groups_max; ii++)
    (*fof_group)[ii].FirstHalo  = NULL;

  if (N_halos<1)
  {
    SID_log("No halos in this file... skipping...", SID_LOG_CLOSE);
    H5Fclose(fin_trees);
    return trees_info;
  }

  // read in the trees
  read_trees_and_catalogs(run_globals, snapshot, fin_trees, *halo, N_halos,
      *fof_group, N_fof_groups, requested_forest_id, N_requested_forests,
      &N_halos_kept, &N_fof_groups_kept);

  // close the tree file
  H5Fclose(fin_trees);

  // if subsampling the trees, then update the trees_info to reflect what we now have
  if(subsample_trees)
  {
    trees_info.n_halos = N_halos_kept;
    trees_info.n_fof_groups = N_fof_groups_kept;
  }

  SID_log("Read %d halos in %d fof_groups.", SID_LOG_COMMENT, trees_info.n_halos, trees_info.n_fof_groups);

  SID_log("...done", SID_LOG_CLOSE);

  return trees_info;
}


