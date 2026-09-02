#include <assert.h>
#include <hdf5_hl.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "dist_func.h"
#include "meraxes.h"

void df_init(distribution_function_t* df, double x_min, double x_max, int bins_per_dex, const char* description)
{
  assert(df != NULL);
  assert(x_max > x_min);
  assert(bins_per_dex > 0);

  df->x_min = x_min;
  df->x_max = x_max;

  // Logarithmic binning: bins per dex (or per magnitude for UVLF/DustyLF)
  df->n_bins = (int)((x_max - x_min) * bins_per_dex);
  df->bin_width = (x_max - x_min) / df->n_bins;

  // Allocate memory for bins and counts
  df->bins = (distribution_bin_t*)malloc(df->n_bins * sizeof(distribution_bin_t));
  df->bin_counts = (double*)calloc(df->n_bins, sizeof(double));
  df->bin_variance = (double*)calloc(df->n_bins, sizeof(double));
  assert(df->bins != NULL);
  assert(df->bin_counts != NULL);
  assert(df->bin_variance != NULL);

  // Store description
  if (description != NULL) {
    strncpy(df->description, description, 255);
    df->description[255] = '\0';
  } else {
    strcpy(df->description, "Distribution Function");
  }

  // Initialize bin centers
  for (int i = 0; i < df->n_bins; i++) {
    df->bins[i].center = x_min + (i + 0.5) * df->bin_width;
    df->bins[i].number_density = 0.0;
    df->bins[i].uncertainty = 0.0;
  }
}

void df_free(distribution_function_t* df)
{
  assert(df != NULL);
  if (df->bins != NULL) {
    free(df->bins);
    df->bins = NULL;
  }
  if (df->bin_counts != NULL) {
    free(df->bin_counts);
    df->bin_counts = NULL;
  }
  if (df->bin_variance != NULL) {
    free(df->bin_variance);
    df->bin_variance = NULL;
  }
}

void df_mpi_reduce(distribution_function_t* df, int mpi_rank, int mpi_size)
{
  assert(df != NULL);

  if (mpi_size > 1) {
    // Sum counts across all processes to rank 0 only
    if (mpi_rank == 0) {
      // Rank 0 uses MPI_IN_PLACE to reduce in-place into df->bin_counts
      MPI_Reduce(MPI_IN_PLACE, df->bin_counts, df->n_bins, MPI_DOUBLE, MPI_SUM, 0, run_globals.mpi_comm);
      MPI_Reduce(MPI_IN_PLACE, df->bin_variance, df->n_bins, MPI_DOUBLE, MPI_SUM, 0, run_globals.mpi_comm);
    } else {
      // Other ranks send their data (recvbuf is ignored for non-root ranks)
      MPI_Reduce(df->bin_counts, NULL, df->n_bins, MPI_DOUBLE, MPI_SUM, 0, run_globals.mpi_comm);
      MPI_Reduce(df->bin_variance, NULL, df->n_bins, MPI_DOUBLE, MPI_SUM, 0, run_globals.mpi_comm);
    }
  }

  if (mpi_rank == 0)
    df_finalize(df);
}

void df_finalize(distribution_function_t* df)
{
  assert(df != NULL);
  assert(df->bins != NULL);

  double total_variance = 0.0;
  for (int i = 0; i < df->n_bins; i++)
    total_variance += df->bin_variance[i];

  const int use_bernoulli = total_variance > 0.0;
  const double normalization = df->volume * df->bin_width;

  for (int i = 0; i < df->n_bins; i++) {
    const double count = df->bin_counts[i];
    df->bins[i].number_density = count / normalization;
    if (use_bernoulli)
      df->bins[i].uncertainty = sqrt(df->bin_variance[i]) / normalization;
    else
      df->bins[i].uncertainty = count > 0.0 ? sqrt(count) / normalization : 0.0;
  }
}

void df_write_hdf5(hid_t file_id,
                   const char* group_name,
                   const distribution_function_t* df,
                   const char* dataset_prefix,
                   const char* units)
{
  assert(df != NULL);
  assert(df->bins != NULL);

  // Open group
  hid_t group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
  if (group_id < 0) {
    mlog_error("Failed to open group %s for writing distribution function", group_name);
    return;
  }

  // Create 2D array: (n_bins, 3) for [center, density, uncertainty]
  hsize_t dims[2] = {(hsize_t)df->n_bins, 3};
  double* data = (double*)malloc(df->n_bins * 3 * sizeof(double));
  assert(data != NULL);

  for (int i = 0; i < df->n_bins; i++) {
    // Column 0: bin center
    data[i * 3 + 0] = df->bins[i].center;
    // Columns 1 and 2 were finalized by df_mpi_reduce/df_finalize.
    data[i * 3 + 1] = df->bins[i].number_density;
    data[i * 3 + 2] = df->bins[i].uncertainty;
  }

  // Write single 2D dataset
  H5LTmake_dataset_double(group_id, dataset_prefix, 2, dims, data);

  // Write metadata as attributes to the dataset
  H5LTset_attribute_int(group_id, dataset_prefix, "n_bins", (int*)&df->n_bins, 1);
  H5LTset_attribute_double(group_id, dataset_prefix, "x_min", (double*)&df->x_min, 1);
  H5LTset_attribute_double(group_id, dataset_prefix, "x_max", (double*)&df->x_max, 1);
  H5LTset_attribute_double(group_id, dataset_prefix, "bin_width", (double*)&df->bin_width, 1);
  H5LTset_attribute_double(group_id, dataset_prefix, "volume", (double*)&df->volume, 1);
  H5LTset_attribute_string(group_id, dataset_prefix, "description", df->description);
  H5LTset_attribute_string(group_id, dataset_prefix, "units", (char*)units);
  
  // Add column labels as attribute
  const char* column_names = "center,density,uncertainty";
  H5LTset_attribute_string(group_id, dataset_prefix, "columns", (char*)column_names);

  H5Gclose(group_id);
  free(data);
}
