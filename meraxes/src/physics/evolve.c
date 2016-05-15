#include <math.h>
#include "meraxes.h"

//! Evolve existing galaxies forward in time
int evolve_galaxies(run_globals_t *run_globals, fof_group_t *fof_group, int snapshot, int NGal, int NFof)
{
  galaxy_t *gal        = NULL;
  halo_t *halo         = NULL;
  int gal_counter      = 0;
  int dead_gals        = 0;
  double infalling_gas = 0;
  double cooling_mass  = 0;
  int NSteps           = run_globals->params.NSteps;
  bool Flag_IRA    = (bool)(run_globals->params.physics.Flag_IRA);

  SID_log("Doing physics...", SID_LOG_OPEN | SID_LOG_TIMER);


  for (int i_fof = 0; i_fof < NFof; i_fof++)
  {
    // First check to see if this FOF group is empty.  If it is then skip it.
    if (fof_group[i_fof].FirstOccupiedHalo == NULL)
      continue;

    infalling_gas = gas_infall(run_globals, &(fof_group[i_fof]), snapshot);

    for (int i_step = 0; i_step < NSteps; i_step++)
    {
      halo = fof_group[i_fof].FirstHalo;
      while (halo != NULL)
      {
        gal = halo->Galaxy;

        while (gal != NULL)
        {
          if (gal->Type == 0)
          {
            cooling_mass = gas_cooling(run_globals, gal);

            add_infall_to_hot(gal, infalling_gas / ((double)NSteps));

            reincorporate_ejected_gas(run_globals, gal);

            cool_gas_onto_galaxy(gal, cooling_mass);
          }

          if (gal->Type < 3)
          {
            if (!Flag_IRA)
            {
              evolve_stellar_pops(run_globals, gal, snapshot);
              delayed_supernova_feedback(run_globals, gal, snapshot);
            }

            insitu_star_formation(run_globals, gal, snapshot);

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
      while (halo != NULL)
      {
        gal = halo->Galaxy;
        while (gal != NULL)
        {
          if (gal->Type == 2)
          {
            // If the merger clock has run out or our target halo has already
            // merged then process a merger event.
            if ((gal->MergTime < 0) || (gal->MergerTarget->Type == 3))
              merge_with_target(run_globals, gal, &dead_gals, snapshot);
            else if (run_globals->params.physics.Flag_BHFeedback)
              previous_merger_driven_BH_growth(run_globals, gal);
          }
          else if (run_globals->params.physics.Flag_BHFeedback)
            previous_merger_driven_BH_growth(run_globals, gal);
          gal = gal->NextGalInHalo;
        }
        halo = halo->NextHaloInFOFGroup;
      }
    }
  }

  if (gal_counter + (run_globals->NGhosts) != NGal)
  {
    SID_log_error("We have not processed the expected number of galaxies...");
    SID_log("gal_counter = %d but NGal = %d", SID_LOG_COMMENT, gal_counter, NGal);
    ABORT(EXIT_FAILURE);
  }

  SID_log("...done", SID_LOG_CLOSE);

  return gal_counter - dead_gals;
}


void passively_evolve_ghost(run_globals_t *run_globals, galaxy_t *gal, int snapshot)
{
  // Passively evolve ghosts.
  // Currently, this just means evolving their stellar pops...

  bool Flag_IRA = (bool)(run_globals->params.physics.Flag_IRA);

  if (!Flag_IRA)
  {
    evolve_stellar_pops(run_globals, gal, snapshot);
    delayed_supernova_feedback(run_globals, gal, snapshot);
  }
}

