#include "emission_lines.h"
#include "core/float_precision_check.h"
#include <math.h>

static double clamp_non_finite(double value)
{
  if (!isfinite(value) || value < 0.0)
    return 0.0;

  return value;
}

void set_OIII_coeffs(double T)
{
  if (T <= 0.0)
    T = 1e4;

  const double beta_q = sqrt(2.0 * M_PI * pow(HBAR_SI, 4) / (BOLTZMANN_SI * pow(ELECTRON_MASS_SI, 3))) * 1e6;
  const double o30 = 0.243 * pow(T / 1e4, 0.120 + 0.031 * log(T / 1e4));
  const double o40 = 0.0321 * pow(T / 1e4, 0.118 + 0.057 * log(T / 1e4));
  const double k03 = (beta_q / sqrt(T)) * o30 * exp(-29169.0 / T);
  const double k04 = (beta_q / sqrt(T)) * o40 * exp(-61207.0 / T);

  run_globals.loiii_params.oxygen_abundance_over_z_sun = (pow(10.0, 8.69) / 1e12) / Z_SUN;
  run_globals.loiii_params.excitation_rate = k03 + k04 * (LOIII_A43 / (LOIII_A43 + LOIII_A41));
  run_globals.loiii_params.branching_ratio = LOIII_A32 / (LOIII_A32 + LOIII_A31);
}

// Adds one star-formation event's intrinsic [O III] luminosity to gal->LOIII_dusty
// (units of 1e40 erg/s). Call from update_reservoirs_from_sf with the event's own
// SFR and the metallicity of the gas that forms the stars, before ColdGas and
// MetalsColdGas are updated. Star formation removes gas at this Z, so it is also the
// cold-gas Z after star formation and before this event's SN feedback.
void compute_LOIII(galaxy_t* gal, double sfr, double metallicity)
{
  if (!(sfr > 0.0))
    return;

  const double sfr_gs = sfr * run_globals.units.UnitMass_in_g / run_globals.units.UnitTime_in_s;

  // The photon rate per superbubble times the number of superbubbles is
  // 4000 SFR / m_p, so the luminosity does not depend on the disc state.
  const double total_ionizing_photon_rate = sfr_gs / PROTONMASS * 4000.0;
  const double oiii_volume_fraction = 0.8;
  // Scaled by 1e-40 so this is stored in units of 1e40 erg/s (matches
  // galaxy_output_t.LOIII_dusty) and stays well within float range.
  const double loiii =
    clamp_non_finite(run_globals.loiii_params.oxygen_abundance_over_z_sun * metallicity *
                     run_globals.loiii_params.excitation_rate * run_globals.loiii_params.branching_ratio *
                     (total_ionizing_photon_rate / ALPHA_HII) * PLANCK * nu32 * oiii_volume_fraction * 1e-40);
  if (loiii <= 0.0)
    return;

  // Stash the intrinsic (non-dust-attenuated) luminosity in LOIII_dusty --
  // there's no separate LOIII field any more. The output/cache preparation
  // code (save.c) reads this as its intrinsic input and overwrites it in
  // place with the real dust-attenuated value once computed.
  gal->LOIII_dusty = check_float_cast((double)(gal->LOIII_dusty) + loiii, FloatField_LOIII_dusty);
}
