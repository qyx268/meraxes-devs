#ifndef FLOAT_PRECISION_CHECK_H
#define FLOAT_PRECISION_CHECK_H

#include <math.h>

#ifdef __cplusplus
extern "C"
{
#endif

// One entry per galaxy_t field that has been converted from double to
// float. Array-valued fields (e.g. NewStars[N_HISTORY_SNAPS]) share a single
// entry across all of their indices.
#define FLOAT_FIELD_LIST(X)                                                                                           \
  X(NewStars)                                                                                                         \
  X(NewStars_II)                                                                                                      \
  X(NewStars_III)                                                                                                     \
  X(NewMetals)                                                                                                        \
  X(inBCFlux)                                                                                                         \
  X(outBCFlux)                                                                                                        \
  X(inBCFluxIII)                                                                                                      \
  X(outBCFluxIII)                                                                                                     \
  X(Mvir)                                                                                                             \
  X(Rvir)                                                                                                             \
  X(Vvir)                                                                                                             \
  X(Vmax)                                                                                                             \
  X(Spin)                                                                                                             \
  X(dt)                                                                                                               \
  X(HotGas)                                                                                                           \
  X(MetalsHotGas)                                                                                                     \
  X(ColdGas)                                                                                                          \
  X(MetalsColdGas)                                                                                                    \
  X(Mcool)                                                                                                            \
  X(StellarMass)                                                                                                      \
  X(GrossStellarMass)                                                                                                 \
  X(Fesc)                                                                                                             \
  X(FescWeightedGSM)                                                                                                  \
  X(FescWeightedSfr)                                                                                                  \
  X(MetalsStellarMass)                                                                                                \
  X(DiskScaleLength)                                                                                                  \
  X(Sfr)                                                                                                              \
  X(LOIII)                                                                                                            \
  X(LOIII_dusty)                                                                                                      \
  X(ionization_param)                                                                                                 \
  X(EjectedGas)                                                                                                       \
  X(MetalsEjectedGas)                                                                                                 \
  X(BlackHoleMass)                                                                                                    \
  X(FescBH)                                                                                                           \
  X(BHemissivity)                                                                                                     \
  X(QuasarLuv)                                                                                                        \
  X(EffectiveBHM)                                                                                                     \
  X(EffectiveBHAR)                                                                                                    \
  X(DutyCycleAGN)                                                                                                     \
  X(BlackHoleAccretedHotMass)                                                                                         \
  X(BlackHoleAccretedColdMass)                                                                                        \
  X(BlackHoleAccretingColdMass)                                                                                       \
  X(BHAccretionOnTime)                                                                                                \
  X(t_resp)                                                                                                          \
  X(SfrIII)                                                                                                           \
  X(StellarMass_II)                                                                                                   \
  X(StellarMass_III)                                                                                                  \
  X(GrossStellarMassIII)                                                                                              \
  X(FescIII)                                                                                                          \
  X(FescIIIWeightedGSM)                                                                                               \
  X(FescIIIWeightedSfr)                                                                                               \
  X(Remnant_Mass)                                                                                                     \
  X(RmetalBubble)                                                                                                     \
  X(PrefactorBubble)                                                                                                  \
  X(TimeBubble)                                                                                                       \
  X(Metal_Probability)                                                                                                \
  X(GalMetal_Probability)                                                                                             \
  X(Metals_IGM)                                                                                                       \
  X(Gas_IGM)                                                                                                          \
  X(Metallicity_IGM)                                                                                                  \
  X(MaxBubble)                                                                                                        \
  X(AveBubble)                                                                                                        \
  X(Prefactor)                                                                                                        \
  X(Times)                                                                                                            \
  X(Radii)                                                                                                            \
  X(mwmsa_num)                                                                                                        \
  X(mwmsa_denom)                                                                                                      \
  X(Rcool)                                                                                                            \
  X(Cos_Inc)                                                                                                          \
  X(MergTime)                                                                                                         \
  X(MergerStartRadius)                                                                                                \
  X(BaryonFracModifier)                                                                                               \
  X(FOFMvirModifier)                                                                                                  \
  X(MvirCrit)                                                                                                         \
  X(MvirCrit_MC)                                                                                                      \
  X(tau_cgm)                                                                                                          \
  X(cumulative_ionization)                                                                                            \
  X(MergerBurstMass)                                                                                                  \
  X(HaloMvir)                                                                                                         \
  X(HaloRvir)                                                                                                         \
  X(HaloVvir)                                                                                                         \
  X(FOFGroupMvir)                                                                                                     \
  X(FOFGroupRvir)                                                                                                     \
  X(FOFGroupVvir)                                                                                                     \
  X(FOFGroupMvirModifier)

#define X_ENUM(name) FloatField_##name,
typedef enum
{
  FLOAT_FIELD_LIST(X_ENUM) FloatField_COUNT
} FloatField;
#undef X_ENUM

#ifdef CHECK_FLOAT_PRECISION

void float_precision_record(FloatField field, double value);
void float_precision_report(void);

static inline float check_float_cast(double value, FloatField field)
{
  float_precision_record(field, value);
  return (float)value;
}

#else

static inline float check_float_cast(double value, FloatField field)
{
  (void)field;
  return (float)value;
}

#define float_precision_report() ((void)0)

#endif // CHECK_FLOAT_PRECISION

#ifdef __cplusplus
}
#endif

#endif // FLOAT_PRECISION_CHECK_H
