#include <math.h>

#include "core/float_precision_check.h"
#include "core/misc_tools.h"
#include "core/modifiers.h"
#include "infall.h"
#include "meraxes.h"
#include "reionization.h"
#include "tree_flags.h"

double gas_infall(fof_group_t* FOFgroup, int snapshot)
{
  halo_t* halo = FOFgroup->FirstHalo;

  galaxy_t* gal;
  galaxy_t* central;
  double total_baryons = 0.;
  double infall_mass = 0.;
  double FOF_Mvir = FOFgroup->Mvir;
  double fb_modifier;

  double total_stellarmass = 0.0;
  double total_hotgas = 0.0;
  double total_coldgas = 0.0;
  double total_ejectedgas = 0.0;
  double total_blackholemass = 0.0;
#if USE_MINI_HALOS 
  double total_remnantmass = 0.0; // This come either from the BHs formed after Pop III stars that fail becoming SN and
                                  // directly collapse or from CCSN. Atm these don't accrete and they don't do anything.
#endif

  // Calculate the total baryon mass in the FOF group
  central = FOFgroup->FirstOccupiedHalo->Galaxy;

  while (halo != NULL) {
    gal = halo->Galaxy;
    while (gal != NULL) {
      total_stellarmass += gal->StellarMass;
      total_hotgas += gal->HotGas;
      total_coldgas += gal->ColdGas;
      total_ejectedgas += gal->EjectedGas;
      total_blackholemass += gal->BlackHoleMass + gal->BlackHoleAccretingColdMass;
#if USE_MINI_HALOS 
      total_remnantmass += gal->Remnant_Mass;
#endif

      if (gal != central) {
        central->HotGas = check_float_cast((double)(central->HotGas) + (gal->HotGas + gal->EjectedGas), FloatField_HotGas);
        central->MetalsHotGas = check_float_cast((double)(central->MetalsHotGas) + (gal->MetalsHotGas + gal->MetalsEjectedGas), FloatField_MetalsHotGas);
        gal->HotGas = check_float_cast((double)(0.0), FloatField_HotGas);
        gal->MetalsHotGas = check_float_cast((double)(0.0), FloatField_MetalsHotGas);
        gal->EjectedGas = check_float_cast((double)(0.0), FloatField_EjectedGas);
        gal->MetalsEjectedGas = check_float_cast((double)(0.0), FloatField_MetalsEjectedGas);
      }

      gal = gal->NextGalInHalo;
    }
    halo = halo->NextHaloInFOFGroup;
  }

  if (check_for_flag(TREE_CASE_BELOW_VIRIAL_THRESHOLD, FOFgroup->FirstHalo->TreeFlags)) {
    // no infall and no hydrstatic hot halo in this case
    central->BaryonFracModifier = check_float_cast((double)(0.0), FloatField_BaryonFracModifier);
    central->EjectedGas = check_float_cast((double)(central->EjectedGas) + (central->HotGas), FloatField_EjectedGas);
    central->MetalsEjectedGas = check_float_cast((double)(central->MetalsEjectedGas) + (central->MetalsHotGas), FloatField_MetalsEjectedGas);
    central->HotGas = check_float_cast((double)(0.0), FloatField_HotGas);
    central->MetalsHotGas = check_float_cast((double)(0.0), FloatField_MetalsHotGas);
    return 0.0;
  }

  total_baryons = total_stellarmass + total_hotgas + total_coldgas + total_ejectedgas + total_blackholemass;
#if USE_MINI_HALOS 
  total_baryons += total_remnantmass;
#endif

  // Calculate the amount of fresh gas required to provide the baryon
  // fraction of this halo.
  fb_modifier = reionization_modifier(central, FOF_Mvir, snapshot);
  if (run_globals.RequestedBaryonFracModifier == 1)
    // N.B. no longer divides by a FOFMvirModifier (fof_group_t no longer has
    // one -- see meraxes.h) -- it was always exactly 1.0 in every real run,
    // since mass-ratio-modifier support has been removed entirely.
    fb_modifier *= interpolate_modifier(run_globals.baryon_frac_modifier,
                                        log10(FOF_Mvir / run_globals.params.Hubble_h) + 10.0);
  infall_mass = fb_modifier * run_globals.params.BaryonFrac * FOF_Mvir - total_baryons;

  // record the infall modifier
  central->BaryonFracModifier = check_float_cast((double)(fb_modifier), FloatField_BaryonFracModifier);

  return infall_mass;
}

void add_infall_to_hot(galaxy_t* central, double infall_mass)
{
#if USE_MINI_HALOS
  bool Flag_Metals = (bool)(run_globals.params.Flag_IncludeMetalEvo);
#endif
  // if we have mass to add then give it to the central
  if (infall_mass > 0) {
    central->HotGas = check_float_cast((double)(central->HotGas) + (infall_mass), FloatField_HotGas);
#if USE_MINI_HALOS
    if (Flag_Metals == true) {
      if (central->Flag_ExtMetEnr == 1) // If the halo is externally enriched, it will accrete polluted gas (metals).
        central->MetalsHotGas = check_float_cast((double)(central->MetalsHotGas) + (infall_mass * central->Metallicity_IGM), FloatField_MetalsHotGas);
    }
#endif
  } else {
    double strip_mass = -infall_mass;
    // otherwise, strip the mass from the ejected
    if (central->EjectedGas > 0) {
      double metallicity = calc_metallicity(central->EjectedGas, central->MetalsEjectedGas);
      central->EjectedGas = check_float_cast((double)(central->EjectedGas) - (strip_mass), FloatField_EjectedGas);
      if (central->EjectedGas < 0) {
        strip_mass = -central->EjectedGas;
        central->EjectedGas = check_float_cast((double)(0.0), FloatField_EjectedGas);
        central->MetalsEjectedGas = check_float_cast((double)(0.0), FloatField_MetalsEjectedGas);
      } else {
        central->MetalsEjectedGas = check_float_cast((double)(central->MetalsEjectedGas) - (metallicity * strip_mass), FloatField_MetalsEjectedGas);
        strip_mass = 0.0;
      }
    }

    // if we still have mass left to remove after exhausting the mass of
    // the ejected component, the remove as much as we can from the hot gas
    if (strip_mass > 0) {
      double metallicity = calc_metallicity(central->HotGas, central->MetalsHotGas);
      central->HotGas = check_float_cast((double)(central->HotGas) - (strip_mass), FloatField_HotGas);
      if (central->HotGas < 0) {
        central->HotGas = check_float_cast((double)(0.0), FloatField_HotGas);
        central->MetalsHotGas = check_float_cast((double)(0.0), FloatField_MetalsHotGas);
      } else
        central->MetalsHotGas = check_float_cast((double)(central->MetalsHotGas) - (metallicity * strip_mass), FloatField_MetalsHotGas);
    }
  }
}
