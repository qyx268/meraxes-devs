#include <assert.h>

#include "float_precision_check.h"
#include "galaxies.h"
#include "magnitudes.h"
#include "meraxes.h"
#include "misc_tools.h"
#include "tree_flags.h"
#include "virial_properties.h"

galaxy_t* new_galaxy(int snapshot, unsigned long halo_ID)
{
  galaxy_t* gal = malloc(sizeof(galaxy_t));

  // Initialise the properties
  gal->ID = (unsigned long)(snapshot * 1e10 + halo_ID);
  gal->Type = -1;
  gal->OldType = -1;
  gal->SnapSkipCounter = 0;
  gal->HaloDescIndex = -1;
  gal->TreeFlags = 0;
  gal->LastIdentSnap = -1;
  gal->Halo = NULL;
  gal->FirstGalInHalo = NULL;
  gal->NextGalInHalo = NULL;
  gal->Next = NULL;
  gal->MergerTarget = NULL;
  gal->Len = 0;
  gal->MaxLen = 0;
  gal->dt = check_float_cast((double)(0.0), FloatField_dt);
  gal->Mvir = check_float_cast((double)(0.0), FloatField_Mvir);
  gal->Rvir = check_float_cast((double)(0.0), FloatField_Rvir);
  gal->Vvir = check_float_cast((double)(0.0), FloatField_Vvir);
  gal->Vmax = check_float_cast((double)(0.0), FloatField_Vmax);
  gal->Spin = check_float_cast((double)(0.0), FloatField_Spin);
  gal->DiskScaleLength = check_float_cast((double)(0.0), FloatField_DiskScaleLength);
  gal->HotGas = check_float_cast((double)(0.0), FloatField_HotGas);
  gal->MetalsHotGas = check_float_cast((double)(0.0), FloatField_MetalsHotGas);
  gal->ColdGas = check_float_cast((double)(0.0), FloatField_ColdGas);
  gal->MetalsColdGas = check_float_cast((double)(0.0), FloatField_MetalsColdGas);
  gal->EjectedGas = check_float_cast((double)(0.0), FloatField_EjectedGas);
  gal->MetalsEjectedGas = check_float_cast((double)(0.0), FloatField_MetalsEjectedGas);
  gal->Mcool = check_float_cast((double)(0.0), FloatField_Mcool);
  gal->Rcool = check_float_cast((double)(0.0), FloatField_Rcool);
  gal->StellarMass = check_float_cast((double)(0.0), FloatField_StellarMass);
  gal->GrossStellarMass = check_float_cast((double)(0.0), FloatField_GrossStellarMass);
  gal->Fesc = check_float_cast((double)(1.0), FloatField_Fesc);
  gal->FescWeightedGSM = check_float_cast((double)(0.0), FloatField_FescWeightedGSM);
  gal->MetalsStellarMass = check_float_cast((double)(0.0), FloatField_MetalsStellarMass);
  gal->LOIII = check_float_cast((double)(0.0), FloatField_LOIII);
  gal->ionization_param = check_float_cast((double)(0.0), FloatField_ionization_param);
  gal->mwmsa_num = check_float_cast((double)(0.0), FloatField_mwmsa_num);
  gal->mwmsa_denom = check_float_cast((double)(0.0), FloatField_mwmsa_denom);
  gal->BlackHoleMass = check_float_cast((double)(run_globals.params.physics.BlackHoleSeed), FloatField_BlackHoleMass);
  gal->FescBH = check_float_cast((double)(1.0), FloatField_FescBH);
  gal->BHemissivity = check_float_cast((double)(0.0), FloatField_BHemissivity);
  gal->QuasarLuv = check_float_cast((double)(0.0), FloatField_QuasarLuv);
  gal->EffectiveBHM = check_float_cast((double)(0.0), FloatField_EffectiveBHM);
  gal->EffectiveBHAR = check_float_cast((double)(0.0), FloatField_EffectiveBHAR);
  gal->DutyCycleAGN = check_float_cast((double)(0.0), FloatField_DutyCycleAGN);
  gal->BlackHoleAccretedHotMass = check_float_cast((double)(0.0), FloatField_BlackHoleAccretedHotMass);
  gal->BlackHoleAccretedColdMass = check_float_cast((double)(0.0), FloatField_BlackHoleAccretedColdMass);
  gal->BlackHoleAccretingColdMass = check_float_cast((double)(0.0), FloatField_BlackHoleAccretingColdMass);
  gal->BHAccretionOnTime = check_float_cast((double)(-1.0), FloatField_BHAccretionOnTime);
  gal->Sfr = check_float_cast((double)(0.0), FloatField_Sfr);
  gal->FescWeightedSfr = check_float_cast((double)(0.0), FloatField_FescWeightedSfr);
  gal->Cos_Inc = check_float_cast((double)(gsl_rng_uniform(run_globals.random_generator)), FloatField_Cos_Inc);
  gal->MergTime = check_float_cast((double)(99999.9), FloatField_MergTime);
  gal->BaryonFracModifier = check_float_cast((double)(1.0), FloatField_BaryonFracModifier);
  gal->FOFMvirModifier = check_float_cast((double)(1.0), FloatField_FOFMvirModifier);
  gal->MvirCrit = check_float_cast((double)(0.0), FloatField_MvirCrit);
  gal->tau_cgm = check_float_cast((double)(0.0), FloatField_tau_cgm);
  gal->cumulative_ionization = check_float_cast((double)(0.0), FloatField_cumulative_ionization);
  gal->MergerBurstMass = check_float_cast((double)(0.0), FloatField_MergerBurstMass);
  gal->MergerStartRadius = check_float_cast((double)(0.0), FloatField_MergerStartRadius);

#if USE_MINI_HALOS
  gal->StellarMass_II = check_float_cast((double)(0.), FloatField_StellarMass_II);
  gal->StellarMass_III = check_float_cast((double)(0.), FloatField_StellarMass_III);
  gal->GrossStellarMassIII = check_float_cast((double)(0.0), FloatField_GrossStellarMassIII);
  gal->SfrIII = check_float_cast((double)(0.0), FloatField_SfrIII);
  gal->FescIII = check_float_cast((double)(1.0), FloatField_FescIII);
  gal->FescIIIWeightedGSM = check_float_cast((double)(0.0), FloatField_FescIIIWeightedGSM);
  gal->FescIIIWeightedSfr = check_float_cast((double)(0.0), FloatField_FescIIIWeightedSfr);
  gal->Remnant_Mass = check_float_cast((double)(0.), FloatField_Remnant_Mass);
  gal->MvirCrit_MC = check_float_cast((double)(0.0), FloatField_MvirCrit_MC);
  gal->Metal_Probability = check_float_cast((double)(0.0), FloatField_Metal_Probability);
  gal->Metals_IGM = check_float_cast((double)(0.0), FloatField_Metals_IGM);
  gal->Gas_IGM = check_float_cast((double)(0.0), FloatField_Gas_IGM);
  gal->Metallicity_IGM = check_float_cast((double)(-50.0), FloatField_Metallicity_IGM);
  gal->RmetalBubble = check_float_cast((double)(0.0), FloatField_RmetalBubble);
  gal->PrefactorBubble = check_float_cast((double)(0.0), FloatField_PrefactorBubble);
  gal->TimeBubble = check_float_cast((double)(0.0), FloatField_TimeBubble);
  gal->AveBubble = check_float_cast((double)(0.), FloatField_AveBubble);
  gal->MaxBubble = check_float_cast((double)(0.), FloatField_MaxBubble);
  gal->Flag_ExtMetEnr = 0;
  gal->GalMetal_Probability = check_float_cast((double)(gsl_rng_uniform(run_globals.random_generator)), FloatField_GalMetal_Probability);

  if (run_globals.params.Flag_IncludeMetalEvo ==
      false) // If you don't have the external metal enrichment all galaxies will start as pristine (Pop.III forming)
    gal->Galaxy_Population = 3;
#endif

  for (int ii = 0; ii < 3; ii++) {
    gal->Pos[ii] = (float)-99999.9;
    gal->Vel[ii] = (float)-99999.9;
  }

  for (int ii = 0; ii < N_HISTORY_SNAPS; ii++) {
    gal->NewStars[ii] = check_float_cast((double)(0.0), FloatField_NewStars);
#if USE_MINI_HALOS
    gal->NewStars_II[ii] = check_float_cast((double)(0.0), FloatField_NewStars_II);
    gal->NewStars_III[ii] = check_float_cast((double)(0.0), FloatField_NewStars_III);
    if (run_globals.params.Flag_IncludeMetalEvo) {
      gal->Prefactor[ii] = check_float_cast((double)(0.0), FloatField_Prefactor);
      gal->Times[ii] = check_float_cast((double)(0.0), FloatField_Times);
      gal->Radii[ii] = check_float_cast((double)(0.0), FloatField_Radii);
    }
#endif
  }

  for (int ii = 0; ii < N_HISTORY_SNAPS; ii++)
    gal->NewMetals[ii] = check_float_cast((double)(0.0), FloatField_NewMetals);

  gal->output_index = -1;
  gal->ghost_flag = false;

#ifdef CALC_MAGS
  init_luminosities(gal);
#endif

  return gal;
}

void copy_halo_props_to_galaxy(halo_t* halo, galaxy_t* gal)
{
  gal->Type = halo->Type;
  gal->Len = halo->Len;
  gal->SnapSkipCounter = halo->SnapOffset;
  gal->HaloDescIndex = halo->DescIndex;
  gal->Mvir = check_float_cast((double)(halo->Mvir), FloatField_Mvir);
  gal->Rvir = check_float_cast((double)(halo->Rvir), FloatField_Rvir);
  gal->Vvir = check_float_cast((double)(halo->Vvir), FloatField_Vvir);
  gal->TreeFlags = halo->TreeFlags;
  gal->Spin = check_float_cast((double)(calculate_spin_param(halo)), FloatField_Spin);
  gal->FOFMvirModifier = check_float_cast((double)(halo->FOFGroup->FOFMvirModifier), FloatField_FOFMvirModifier);

  double sqrt_2 = 1.414213562;
  if (gal->Type == 0) {
    gal->Vmax = check_float_cast((double)(halo->Vmax), FloatField_Vmax);
    gal->DiskScaleLength = check_float_cast((double)(gal->Spin * gal->Rvir / sqrt_2), FloatField_DiskScaleLength);
  } else {
    if (!run_globals.params.physics.Flag_FixVmaxOnInfall)
      gal->Vmax = check_float_cast((double)(halo->Vmax), FloatField_Vmax);
    if (!run_globals.params.physics.Flag_FixDiskRadiusOnInfall)
      gal->DiskScaleLength = check_float_cast((double)(gal->Spin * gal->Rvir / sqrt_2), FloatField_DiskScaleLength);
  }

  for (int ii = 0; ii < 3; ii++) {
    gal->Pos[ii] = halo->Pos[ii];
    gal->Vel[ii] = halo->Vel[ii];
  }

  // record the maximum Len value if necessary
  if (halo->Len > gal->MaxLen)
    gal->MaxLen = halo->Len;
}

void reset_galaxy_properties(galaxy_t* gal, int snapshot)
{
  // Here we reset any galaxy properties which are calculated on a snapshot by
  // snapshot basis.
  gal->Sfr = check_float_cast((double)(0.0), FloatField_Sfr);
  gal->LOIII = check_float_cast((double)(0.0), FloatField_LOIII);
  gal->ionization_param = check_float_cast((double)(0.0), FloatField_ionization_param);
  gal->FescWeightedSfr = check_float_cast((double)(0.0), FloatField_FescWeightedSfr);
  gal->Mcool = check_float_cast((double)(0.0), FloatField_Mcool);
  gal->Rcool = check_float_cast((double)(0.0), FloatField_Rcool);
  gal->MvirCrit = check_float_cast((double)(0.0), FloatField_MvirCrit);
  gal->tau_cgm = check_float_cast((double)(0.0), FloatField_tau_cgm);
  gal->cumulative_ionization = check_float_cast((double)(0.0), FloatField_cumulative_ionization);
  gal->BHemissivity = check_float_cast((double)(0.0), FloatField_BHemissivity);
  gal->QuasarLuv = check_float_cast((double)(0.0), FloatField_QuasarLuv);
  gal->BaryonFracModifier = check_float_cast((double)(1.0), FloatField_BaryonFracModifier);
  gal->FOFMvirModifier = check_float_cast((double)(1.0), FloatField_FOFMvirModifier);
  gal->EffectiveBHAR = check_float_cast((double)(0.0), FloatField_EffectiveBHAR);
  gal->DutyCycleAGN = check_float_cast((double)(0.0), FloatField_DutyCycleAGN);
  gal->BlackHoleAccretedHotMass = check_float_cast((double)(0.0), FloatField_BlackHoleAccretedHotMass);
  gal->BlackHoleAccretedColdMass = check_float_cast((double)(0.0), FloatField_BlackHoleAccretedColdMass);
  gal->t_resp = check_float_cast((double)(1e30), FloatField_t_resp);
#if USE_MINI_HALOS
  gal->SfrIII = check_float_cast((double)(0.0), FloatField_SfrIII);
  gal->FescIIIWeightedSfr = check_float_cast((double)(0.0), FloatField_FescIIIWeightedSfr);
  gal->MvirCrit_MC = check_float_cast((double)(0.0), FloatField_MvirCrit_MC);
#endif

  // Update the stellar mass weighted mean age values.  This only needs to be
  // done for snapshots shich are passing out of what we are able to track
  // with N_HISTORY_SNAPS.
  assert(snapshot > 0);
  if (snapshot >= N_HISTORY_SNAPS) {
    gal->mwmsa_denom = check_float_cast((double)(gal->mwmsa_denom) + (gal->NewStars[N_HISTORY_SNAPS - 1]), FloatField_mwmsa_denom);
    gal->mwmsa_num = check_float_cast((double)(gal->mwmsa_num) + (gal->NewStars[N_HISTORY_SNAPS - 1] * run_globals.LTTime[snapshot - N_HISTORY_SNAPS]), FloatField_mwmsa_num);
  }

  // roll over the baryonic history arrays
  for (int ii = N_HISTORY_SNAPS - 1; ii > 0; ii--) {
    gal->NewStars[ii] = check_float_cast((double)(gal->NewStars[ii - 1]), FloatField_NewStars);
#if USE_MINI_HALOS
    gal->NewStars_II[ii] = check_float_cast((double)(gal->NewStars_II[ii - 1]), FloatField_NewStars_II);
    gal->NewStars_III[ii] = check_float_cast((double)(gal->NewStars_III[ii - 1]), FloatField_NewStars_III);
#endif
  }

  for (int ii = N_HISTORY_SNAPS - 1; ii > 0; ii--)
    gal->NewMetals[ii] = check_float_cast((double)(gal->NewMetals[ii - 1]), FloatField_NewMetals);

  gal->NewStars[0] = check_float_cast((double)(0.0), FloatField_NewStars);
#if USE_MINI_HALOS
  gal->NewStars_II[0] = check_float_cast((double)(0.0), FloatField_NewStars_II);
  gal->NewStars_III[0] = check_float_cast((double)(0.0), FloatField_NewStars_III);
#endif
  gal->NewMetals[0] = check_float_cast((double)(0.0), FloatField_NewMetals);
}

static void push_galaxy_to_halo(galaxy_t* gal, halo_t* halo)
{
  if (halo->Galaxy == NULL)
    halo->Galaxy = gal;
  else {
    // Walk the galaxy list for the halo to find the end and then link the
    // new galaxy
    galaxy_t* cur_gal = halo->Galaxy;
    galaxy_t* prev_gal = cur_gal;
    while (cur_gal != NULL) {
      prev_gal = cur_gal;
      cur_gal->Halo = halo;
      cur_gal = cur_gal->NextGalInHalo;
    }
    prev_gal->NextGalInHalo = gal;
  }

  // Loop through the new galaxy, and all galaxies attached to it, and set
  // the first galaxy in halo and halo pointer
  gal = halo->Galaxy;
  while (gal != NULL) {
    gal->Halo = halo;
    gal->FirstGalInHalo = halo->Galaxy;
    gal = gal->NextGalInHalo;
  }
}

void connect_galaxy_and_halo(galaxy_t* gal, halo_t* halo, int* merger_counter)
{

  if (halo->Galaxy == NULL)
    push_galaxy_to_halo(gal, halo);
  else {
    // There is already a galaxy been assigned to this halo.  That means we
    // have a merger. Now we need to work out which galaxy is merging into
    // which.

    assert(merger_counter != NULL);
    (*merger_counter)++;

    galaxy_t* parent = NULL;
    galaxy_t* infaller = NULL;

    switch (run_globals.params.TreesID) {
      case GBPTREES_TREES:
        // For gbpTrees, we have the merger flags to give us guidance.  Let's use them...
        // TODO: Make sure I don't need to turn off the merger flag...
        parent = check_for_flag(TREE_CASE_MERGER, gal->TreeFlags) ? halo->Galaxy : gal;
        infaller = halo->Galaxy == parent ? gal : halo->Galaxy;
        assert(check_for_flag(TREE_CASE_MERGER, infaller->TreeFlags));
        break;

      case VELOCIRAPTOR_TREES:
      case VELOCIRAPTOR_TREES_AUG:
      case VELOCIRAPTOR_TREES_TRUNCATED:
        // For VELOCIraptor we have some guidance in the form of the progenitor indices.

        if (check_for_flag(TREE_CASE_NO_PROGENITORS, halo->TreeFlags)) {
          // The host halo has been marked as having no progenitors.
          // Since we have a merger though, their are clearly halos which
          // think this is their descendant.  In this case, none of the
          // halo progenitors are deemed to be good enough matches and so
          // we can't use the pointers to select the main progenitor.
          // Let's use the mass instead in this case.
          //
          // N.B. We haven't yet copied any halo properties into the
          // galaxies and so the galaxy.Mvir values still correspond to
          // the previous snapshot.
          parent = gal->Mvir > halo->Galaxy->Mvir ? gal : halo->Galaxy;
        } else {
          // Here we can use the pointers.
          parent = gal->HaloDescIndex == halo->ProgIndex ? gal : halo->Galaxy;
        }

        infaller = halo->Galaxy == parent ? gal : halo->Galaxy;
        break;

      default:
        mlog_error("Unrecognised input trees identifier (TreesID).");
        break;
    }

    infaller->Type = 2;
    // Make sure the halo is pointing to the right galaxy
    if (parent != halo->Galaxy)
      halo->Galaxy = parent;

    // Add the incoming galaxy to the end of the halo's linked list
    push_galaxy_to_halo(infaller, halo);
  }
}

void create_new_galaxy(int snapshot, halo_t* halo, int* NGal, int* new_gal_counter, int* merger_counter)
{
  galaxy_t* gal;

  gal = new_galaxy(snapshot, halo->ID);
  gal->Halo = halo;

  if (snapshot > 0)
    gal->LastIdentSnap = snapshot - 1;
  else
    gal->LastIdentSnap = snapshot;

  connect_galaxy_and_halo(gal, halo, merger_counter);

  if (run_globals.LastGal != NULL)
    run_globals.LastGal->Next = gal;
  else
    run_globals.FirstGal = gal;

  run_globals.LastGal = gal;
  gal->dt = check_float_cast((double)(run_globals.LTTime[gal->LastIdentSnap] - run_globals.LTTime[snapshot]), FloatField_dt);
  *NGal = *NGal + 1;
  *new_gal_counter = *new_gal_counter + 1;
}

void kill_galaxy(galaxy_t* gal, galaxy_t* prev_gal, int* NGal, int* kill_counter)
{
  galaxy_t* cur_gal;

  // Remove it from the global linked list
  if (prev_gal != NULL)
    prev_gal->Next = gal->Next;
  else
    run_globals.FirstGal = gal->Next;

  cur_gal = gal->FirstGalInHalo;

  if (cur_gal != gal) {
    // If it is a type 2 then also remove it from the linked list of galaxies in its halo
    while ((cur_gal->NextGalInHalo != NULL) && (cur_gal->NextGalInHalo != gal))
      cur_gal = cur_gal->NextGalInHalo;
    cur_gal->NextGalInHalo = gal->NextGalInHalo;
  } else {
    // If it is a type 0 or 1 (i.e. first galaxy in it's halo) and there are
    // other galaxies in this halo, reset the FirstGalInHalo pointer so that
    // the satellites can be killed later
    cur_gal = gal->NextGalInHalo;
    while (cur_gal != NULL) {
      cur_gal->FirstGalInHalo = gal->NextGalInHalo;
      cur_gal = cur_gal->NextGalInHalo;
    }
  }

  // Finally deallocated the galaxy and decrement any necessary counters
  free(gal);
  *NGal = *NGal - 1;
  *kill_counter = *kill_counter + 1;
}
