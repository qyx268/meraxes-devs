#include "evolve.h"
#include "blackhole_feedback.h"
#include "cooling.h"
#include "core/stellar_feedback.h"
#if USE_MINI_HALOS
#include "core/PopIII.h"
#include "core/misc_tools.h"
#include "core/virial_properties.h"
#endif
#include "infall.h"
#include "meraxes.h"
#include "mergers.h"
#include "reincorporation.h"
#include "star_formation.h"
#include "supernova_feedback.h"
#include <math.h>

//! Evolve existing galaxies forward in time
#if USE_MINI_HALOS
int evolve_galaxies(fof_group_t* fof_group,
                    int snapshot,
                    int NGal,
                    int NFof,
                    int* gal_counter_Pop3,
                    int* gal_counter_Pop2,
                    int* gal_counter_enriched)
#else
int evolve_galaxies(fof_group_t* fof_group, int snapshot, int NGal, int NFof)
#endif
{
  galaxy_t* gal = NULL;
  halo_t* halo = NULL;
  int gal_counter = 0;
  int dead_gals = 0;
  double infalling_gas = 0;
  double cooling_mass = 0;
  int NSteps = run_globals.params.NSteps;
  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);
#if USE_MINI_HALOS
  double DiskMetallicity; // Need this to compute the internal enrichment in a more accurate way
  bool Flag_Metals = (bool)(run_globals.params.Flag_IncludeMetalEvo);
#endif

  mlog("Doing physics...", MLOG_OPEN | MLOG_TIMERSTART);
  // pre-calculate feedback tables for each lookback snapshot
  compute_stellar_feedback_tables(snapshot);

  for (int i_fof = 0; i_fof < NFof; i_fof++) {
    // First check to see if this FOF group is empty.  If it is then skip it.
    if (fof_group[i_fof].FirstOccupiedHalo == NULL)
      continue;

    infalling_gas = gas_infall(&(fof_group[i_fof]), snapshot);

    for (int i_step = 0; i_step < NSteps; i_step++) {
      halo = fof_group[i_fof].FirstHalo;
      while (halo != NULL) {
        gal = halo->Galaxy;

        while (gal != NULL) {

#if USE_MINI_HALOS
          if (Flag_Metals ==
              true) { // Assign to newly formed galaxies metallicity of their cell according to a certain probability
            if ((gal->Type == 0) &&
                (gal->Flag_ExtMetEnr ==
                 0)) { // In order to be consistent with the rest of Meraxes do this only for the central galaxies!
              if ((gal->GalMetal_Probability <= gal->Metal_Probability) ||
                  (gal->GrossStellarMass + gal->GrossStellarMassIII) > 1e-10) {
                gal->Flag_ExtMetEnr = 1; // Just update the flag. Here what I am saying is that a galaxy that already
                                         // experienced SN events will surely be inside a metal bubble!

                *gal_counter_enriched = *gal_counter_enriched + 1;
                if ((gal->Metallicity_IGM / 0.01) > run_globals.params.physics.ZCrit) {
                  *gal_counter_Pop2 = *gal_counter_Pop2 + 1;
                  gal->Galaxy_Population = 2;
                } else
                  gal->Galaxy_Population = 3; // Enriched but not enough
              }

              else {
                gal->Galaxy_Population = 3;
                gal->Flag_ExtMetEnr = 0;
                *gal_counter_Pop3 = *gal_counter_Pop3 + 1;
              }
            }
          }
#endif

          if (gal->Type == 0) {
            cooling_mass = gas_cooling(gal);

            add_infall_to_hot(
              gal, infalling_gas / ((double)NSteps)); // This function is now updated! If the gal is externally
                                                      // enriched, we will add MetalHotGas according to IGM metallicity!

            reincorporate_ejected_gas(gal);

            cool_gas_onto_galaxy(gal, cooling_mass);
          }

          if (gal->Type < 3) {
            if (!Flag_IRA)
              delayed_supernova_feedback(gal, snapshot);

            if (gal->BlackHoleAccretingColdMass > 0)
              previous_merger_driven_BH_growth(gal);

#if USE_MINI_HALOS
            DiskMetallicity = calc_metallicity(
              gal->ColdGas, gal->MetalsColdGas); // A more accurate way to account for the internal enrichment!
            if ((DiskMetallicity / 0.01) > run_globals.params.physics.ZCrit)
              gal->Galaxy_Population = 2;
            else
              gal->Galaxy_Population = 3;
#endif

            insitu_star_formation(gal, snapshot);

#if USE_MINI_HALOS
            if ((Flag_Metals == true) && (gal->Type < 3)) { // For gal->Type > 0 you are just letting the bubble grow
              calc_metal_bubble(gal, snapshot);
            }
#endif
            // If this is a type 2 then decrement the merger clock
            if (gal->Type == 2)
              gal->MergTime -= gal->dt;
          }

          if (i_step == NSteps - 1)
            gal_counter++;

          gal = gal->NextGalInHalo;
        }

        halo = halo->NextHaloInFOFGroup;
      }

      // Check for mergers
      halo = fof_group[i_fof].FirstHalo;
      while (halo != NULL) {
        gal = halo->Galaxy;
        while (gal != NULL) {
          if (gal->Type == 2)
            // If the merger clock has run out or our target halo has already
            // merged then process a merger event.
            if ((gal->MergTime < 0) || (gal->MergerTarget->Type == 3))
              merge_with_target(gal, &dead_gals, snapshot);

          gal = gal->NextGalInHalo;
        }
        halo = halo->NextHaloInFOFGroup;
      }
    }
  }

  if (gal_counter + (run_globals.NGhosts) != NGal) {
    mlog_error("We have not processed the expected number of galaxies...");
    mlog("gal_counter = %d but NGal = %d", MLOG_MESG, gal_counter, NGal);
    ABORT(EXIT_FAILURE);
  }

  mlog("...done", MLOG_CLOSE | MLOG_TIMERSTOP);

  return gal_counter - dead_gals;
}

void passively_evolve_ghost(galaxy_t* gal, int snapshot)
{
  // Passively evolve ghosts.
  // Currently, this just means evolving their stellar pops...

  bool Flag_IRA = (bool)(run_globals.params.physics.Flag_IRA);
#if USE_MINI_HALOS
  bool Flag_Metals = (bool)(run_globals.params.Flag_IncludeMetalEvo);
#endif

  if (!Flag_IRA)
    delayed_supernova_feedback(gal, snapshot);
}
