#ifndef SAVE_H
#define SAVE_H

#include "meraxes.h"

typedef struct galaxy_output_t
{
#ifdef CALC_MAGS
  float Mags[MAGS_N_BANDS];
  float DustyMags[MAGS_N_BANDS];
#if USE_MINI_HALOS
  float MagsIII[MAGS_N_BANDS];
#endif
#endif

  int Type;
  int CentralGal;
  int GhostFlag;
  int Len;

  float Pos[3];
  float Spin;
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float FOFMvir;

  // baryonic reservoirs
  float HotGas;
  float MetalsHotGas;
  float ColdGas;
  float MetalsColdGas;
  float Mcool;
  float DiskScaleLength;
  float StellarMass;
  float GrossStellarMass;
  float Fesc;
  float FescWeightedGSM;
  float FescWeightedSfr;
  float MetalsStellarMass;
  float Sfr;
  float LOIII;
#ifdef CALC_MAGS
  float LOIII_dusty;
#endif
  float EjectedGas;
  float MetalsEjectedGas;
  float BlackHoleMass;
  float BHemissivity;
  float QuasarMag;
  float EffectiveBHM;
  float DutyCycleAGN;

#if USE_MINI_HALOS
  int Galaxy_Population; // Disentangles Pop III/Pop II; not needed when Mini_halos is off, since all galaxies are Pop II
  float GrossStellarMassIII;
  float FescIII;
  float FescIIIWeightedGSM;
  float FescIIIWeightedSfr;

  float RmetalBubble;
  int Flag_ExtMetEnr;
  float Metal_Probability;
  float GalMetal_Probability;
  float StellarMass_II;
  float StellarMass_III;
  float Remnant_Mass;
  float MvirCrit_MC;
#endif

  // misc
  float Rcool;
  float MergTime;
  float BaryonFracModifier;
  float MvirCrit;
  float tau_cgm;
  float dt;

  // baryonic histories
  float MWMSA; // Mass weighted mean stellar age
  float NewStars[N_HISTORY_SNAPS];
#if USE_MINI_HALOS
  float NewStars_II[N_HISTORY_SNAPS];
  float NewStars_III[N_HISTORY_SNAPS];
#endif
} galaxy_output_t;

#ifdef __cplusplus
extern "C"
{
#endif

  void prepare_galaxy_for_output(struct galaxy_t gal, galaxy_output_t* galout, int i_snap);
  void calc_hdf5_props(void);
  void prep_hdf5_file(void);
  void create_master_file(void);
  void prepare_luminosity_function_cache(int snapshot);
  void write_snapshot(int n_write, int i_out, int* last_n_write);

#ifdef __cplusplus
}
#endif

#endif
