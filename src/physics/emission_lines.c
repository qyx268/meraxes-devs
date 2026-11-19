#include "emission_lines.h"
#include "core/float_precision_check.h"
#include <math.h>

static double clamp_non_finite(double value)
{
  if (!isfinite(value) || value < 0.0)
    return 0.0;

  return value;
}

// Metallicity range [Z_sun] of the Yang et al. (2023) electron-temperature fit.
// Outside it the temperature stays at its value at the nearest end.
#define OIII_TEMPERATURE_FIT_Z_MIN 0.016
#define OIII_TEMPERATURE_FIT_Z_MAX 1.6

void set_OIII_coeffs(void)
{
  // sqrt(2 pi hbar^4 / (k m_e^3)) = 8.63e-6 cm^3 s^-1 K^1/2
  run_globals.loiii_params.excitation_prefactor =
    sqrt(2.0 * M_PI * pow(HBAR_SI, 4) / (BOLTZMANN_SI * pow(ELECTRON_MASS_SI, 3))) * 1e6;
  run_globals.loiii_params.oxygen_abundance_over_z_sun = (pow(10.0, 8.69) / 1e12) / Z_SUN;
  run_globals.loiii_params.branching_ratio = LOIII_A32 / (LOIII_A32 + LOIII_A31);
}

// Electron temperature [K] of HII regions in gas of this metallicity, from the
// Yang et al. (2023) fit T4 = 0.824 x^2 + 0.101 x + 1.08 with x = log10(Z / Z_sun).
double oiii_electron_temperature(double metallicity)
{
  // fmax also takes a zero, negative or NaN metallicity to the low end of the fit.
  const double z_ratio = fmin(fmax(metallicity / Z_SUN, OIII_TEMPERATURE_FIT_Z_MIN), OIII_TEMPERATURE_FIT_Z_MAX);
  const double x = log10(z_ratio);
  return 1e4 * (0.824 * x * x + 0.101 * x + 1.08);
}

// Case-B hydrogen recombination coefficient [cm^3/s] (Osterbrock 1989).
double case_b_recombination_coefficient(double temperature)
{
  return 2.59e-13 * pow(temperature / 1e4, -0.7);
}

// Rate coefficient [cm^3/s] for electrons exciting O++ into 1D2, directly from the
// ground term and through 1S0 -> 1D2 cascades; collision strengths from Draine (2011).
double oiii_excitation_rate(double temperature)
{
  const double t4 = temperature / 1e4;
  const double prefactor = run_globals.loiii_params.excitation_prefactor / sqrt(temperature);
  const double k03 = prefactor * 0.243 * pow(t4, 0.120 + 0.031 * log(t4)) * exp(-29169.0 / temperature);
  const double k04 = prefactor * 0.0321 * pow(t4, 0.118 + 0.057 * log(t4)) * exp(-62134.0 / temperature);
  return k03 + k04 * (LOIII_A43 / (LOIII_A43 + LOIII_A41));
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

  // Ionizing photons per stellar baryon as in reionization (Pop. II). The photon
  // rate per superbubble times the number of superbubbles is this total, so the
  // luminosity does not depend on the disc state.
  const double photon_rate = sfr_gs / PROTONMASS * run_globals.params.physics.ReionNionPhotPerBary;
  // The HII regions' temperature, and with it the excitation and recombination
  // rates, follows the metallicity of the gas the stars form from.
  const double temperature = oiii_electron_temperature(metallicity);
  const double alpha_b = case_b_recombination_coefficient(temperature);
  const double oiii_volume_fraction = 0.8;
  // Scaled by 1e-40 so this is stored in units of 1e40 erg/s (matches
  // galaxy_output_t.LOIII_dusty) and stays well within float range.
  const double loiii = clamp_non_finite(run_globals.loiii_params.oxygen_abundance_over_z_sun * metallicity *
                                        oiii_excitation_rate(temperature) * run_globals.loiii_params.branching_ratio *
                                        (photon_rate / alpha_b) * PLANCK * nu32 * oiii_volume_fraction * 1e-40);
  if (loiii <= 0.0)
    return;

  // Stash the intrinsic (non-dust-attenuated) luminosity in LOIII_dusty --
  // there's no separate LOIII field any more. The output/cache preparation
  // code (save.c) reads this as its intrinsic input and overwrites it in
  // place with the real dust-attenuated value once computed.
  gal->LOIII_dusty = check_float_cast((double)(gal->LOIII_dusty) + loiii, FloatField_LOIII_dusty);
}
