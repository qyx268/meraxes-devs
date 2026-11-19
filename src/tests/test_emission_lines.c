#define _MAIN
#include <criterion/criterion.h>
#include <math.h>
#include <string.h>

#include "core/init.h"
#include "meraxes.h"
#include "physics/emission_lines.h"

static void reset_run_globals(void)
{
  memset(&run_globals, 0, sizeof(run_globals));
}

static void configure_default_units_and_params(void)
{
  run_globals.params.Hubble_h = 0.6774;

  // Typical cosmological code units used by existing tests.
  run_globals.units.UnitLength_in_cm = 3.085678e24;   // 1 Mpc/h
  run_globals.units.UnitMass_in_g = 1.989e43;         // 1e10 Msun/h
  run_globals.units.UnitVelocity_in_cm_per_s = 1.0e5; // 1 km/s

  set_units();
}

static galaxy_t make_valid_test_galaxy(void)
{
  galaxy_t gal;
  memset(&gal, 0, sizeof(gal));

  gal.Mvir = 10.0;
  gal.Rvir = 0.08;
  gal.Vmax = 150.0;
  gal.Spin = 0.04;
  gal.dt = 0.01;

  gal.StellarMass = 0.5;
  gal.ColdGas = 0.7;
  gal.MetalsColdGas = 0.01;
  gal.Mcool = 0.2;
  gal.Sfr = 0.5;
  gal.DiskScaleLength = 0.003;

  return gal;
}

TestSuite(emission_lines);

Test(emission_lines, set_OIII_coeffs_generates_finite_coefficients)
{
  reset_run_globals();
  set_OIII_coeffs(1e4);

  cr_expect(isfinite(run_globals.loiii_params.oxygen_abundance_over_z_sun));
  cr_expect(isfinite(run_globals.loiii_params.excitation_rate));
  cr_expect(isfinite(run_globals.loiii_params.branching_ratio));

  cr_expect_gt(run_globals.loiii_params.oxygen_abundance_over_z_sun, 0.0);
  cr_expect_gt(run_globals.loiii_params.excitation_rate, 0.0);
  cr_expect_gt(run_globals.loiii_params.branching_ratio, 0.0);
}

Test(emission_lines, set_OIII_coeffs_non_positive_temperature_falls_back_to_1e4)
{
  reset_run_globals();

  set_OIII_coeffs(1e4);
  const loiii_params_t ref = run_globals.loiii_params;

  set_OIII_coeffs(0.0);

  cr_expect_float_eq(run_globals.loiii_params.oxygen_abundance_over_z_sun, ref.oxygen_abundance_over_z_sun, 1e-15);
  cr_expect_float_eq(run_globals.loiii_params.excitation_rate, ref.excitation_rate, 1e-15);
  cr_expect_float_eq(run_globals.loiii_params.branching_ratio, ref.branching_ratio, 1e-15);
}

static void setup_loiii(void)
{
  reset_run_globals();
  configure_default_units_and_params();
  set_OIII_coeffs(1e4);
}

static const double event_sfr = 0.5;
static const double event_metallicity = 0.01 / 0.7;

// Intrinsic luminosity [1e40 erg/s] of one event on a galaxy without earlier [O III]
// this snapshot.
static double single_event(double sfr, double metallicity)
{
  galaxy_t gal = make_valid_test_galaxy();
  compute_LOIII(&gal, sfr, metallicity);
  return gal.LOIII_dusty;
}

Test(emission_lines, compute_LOIII_valid_event_gives_finite_positive_luminosity)
{
  setup_loiii();
  const double loiii = single_event(event_sfr, event_metallicity);
  cr_expect(isfinite(loiii));
  cr_expect_gt(loiii, 0.0);
}

// Each star-formation event adds its own luminosity, from its own SFR.
Test(emission_lines, compute_LOIII_adds_events)
{
  setup_loiii();
  const double first = single_event(event_sfr, event_metallicity);
  const double second = single_event(3.0 * event_sfr, 0.5 * event_metallicity);

  galaxy_t gal = make_valid_test_galaxy();
  compute_LOIII(&gal, event_sfr, event_metallicity);
  compute_LOIII(&gal, 3.0 * event_sfr, 0.5 * event_metallicity);
  cr_expect_float_eq(gal.LOIII_dusty, first + second, 1e-6 * (first + second));
}

Test(emission_lines, compute_LOIII_is_linear_in_metallicity_and_sfr)
{
  setup_loiii();
  const double reference = single_event(event_sfr, event_metallicity);
  cr_expect_float_eq(single_event(event_sfr, 2.0 * event_metallicity), 2.0 * reference, 1e-6 * reference);
  cr_expect_float_eq(single_event(3.0 * event_sfr, event_metallicity), 3.0 * reference, 1e-6 * reference);
}

// The photon rate per superbubble times the number of superbubbles cancels, so the
// disc state does not enter, even when it is invalid.
Test(emission_lines, compute_LOIII_does_not_depend_on_disc_state)
{
  setup_loiii();
  const double reference = single_event(event_sfr, event_metallicity);

  galaxy_t gal = make_valid_test_galaxy();
  gal.Mvir = 0.0;
  gal.ColdGas = 0.0;
  gal.MetalsColdGas = 0.0;
  gal.DiskScaleLength = 0.0;
  compute_LOIII(&gal, event_sfr, event_metallicity);
  cr_expect_float_eq(gal.LOIII_dusty, reference, 1e-6 * reference);
}

Test(emission_lines, compute_LOIII_without_sfr_or_metals_adds_nothing)
{
  setup_loiii();
  galaxy_t gal = make_valid_test_galaxy();
  gal.LOIII_dusty = 123.0f;

  compute_LOIII(&gal, 0.0, event_metallicity);
  compute_LOIII(&gal, event_sfr, 0.0);
  cr_expect_eq(gal.LOIII_dusty, 123.0f);
}
