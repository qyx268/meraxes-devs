#ifdef CHECK_FLOAT_PRECISION

#include <float.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>

#include "float_precision_check.h"
#include "meraxes.h"

#define NEAR_OVERFLOW_FRACTION 0.01 // "close to overflow" := within 2 orders of magnitude of FLT_MAX

typedef struct
{
  double min;
  double max;
  long n_calls;
  long n_overflow;
  long n_near_overflow;
} float_field_stats_t;

#define X_NAME(name) #name,
static const char* field_names[] = { FLOAT_FIELD_LIST(X_NAME) };
#undef X_NAME

static float_field_stats_t stats[FloatField_COUNT];
static bool initialised = false;

static void ensure_initialised(void)
{
  if (initialised)
    return;

  for (int i = 0; i < FloatField_COUNT; i++) {
    stats[i].min = DBL_MAX;
    stats[i].max = -DBL_MAX;
    stats[i].n_calls = 0;
    stats[i].n_overflow = 0;
    stats[i].n_near_overflow = 0;
  }
  initialised = true;
}

void float_precision_record(FloatField field, double value)
{
  ensure_initialised();

  if (!isfinite(value))
    return;

  float_field_stats_t* s = &stats[field];
  s->n_calls++;
  if (value < s->min)
    s->min = value;
  if (value > s->max)
    s->max = value;

  double magnitude = fabs(value);
  if (magnitude > (double)FLT_MAX)
    s->n_overflow++;
  else if (magnitude > NEAR_OVERFLOW_FRACTION * (double)FLT_MAX)
    s->n_near_overflow++;
}

void float_precision_report(void)
{
  ensure_initialised();

  int mpi_rank = run_globals.mpi_rank;
  MPI_Comm comm = run_globals.mpi_comm;

  double global_min[FloatField_COUNT];
  double global_max[FloatField_COUNT];
  long global_n_calls[FloatField_COUNT];
  long global_n_overflow[FloatField_COUNT];
  long global_n_near_overflow[FloatField_COUNT];

  double local_min[FloatField_COUNT];
  double local_max[FloatField_COUNT];
  long local_n_calls[FloatField_COUNT];
  long local_n_overflow[FloatField_COUNT];
  long local_n_near_overflow[FloatField_COUNT];

  for (int i = 0; i < FloatField_COUNT; i++) {
    local_min[i] = stats[i].min;
    local_max[i] = stats[i].max;
    local_n_calls[i] = stats[i].n_calls;
    local_n_overflow[i] = stats[i].n_overflow;
    local_n_near_overflow[i] = stats[i].n_near_overflow;
  }

  MPI_Reduce(local_min, global_min, FloatField_COUNT, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(local_max, global_max, FloatField_COUNT, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(local_n_calls, global_n_calls, FloatField_COUNT, MPI_LONG, MPI_SUM, 0, comm);
  MPI_Reduce(local_n_overflow, global_n_overflow, FloatField_COUNT, MPI_LONG, MPI_SUM, 0, comm);
  MPI_Reduce(local_n_near_overflow, global_n_near_overflow, FloatField_COUNT, MPI_LONG, MPI_SUM, 0, comm);

  if (mpi_rank != 0)
    return;

  mlog("=== CHECK_FLOAT_PRECISION report (galaxy_t double->float conversion) ===", MLOG_MESG);
  mlog("field                          n_calls     min             max             n_overflow  n_near_overflow",
       MLOG_MESG);
  for (int i = 0; i < FloatField_COUNT; i++) {
    if (global_n_calls[i] == 0)
      continue;
    mlog("%-30s %-11ld %-15.6e %-15.6e %-11ld %-11ld",
         MLOG_MESG,
         field_names[i],
         global_n_calls[i],
         global_min[i],
         global_max[i],
         global_n_overflow[i],
         global_n_near_overflow[i]);
  }
  mlog("=== end CHECK_FLOAT_PRECISION report ===", MLOG_MESG);
}

#endif // CHECK_FLOAT_PRECISION
