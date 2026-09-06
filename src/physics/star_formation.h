#ifndef STAR_FORMATION_H
#define STAR_FORMATION_H

#include "meraxes.h"

typedef enum SFtype
{
  INSITU,
  MERGER
} SFtype;

#ifdef __cplusplus
extern "C"
{
#endif

  void update_reservoirs_from_sf(struct galaxy_t* gal, double new_stars, int snapshot, SFtype type);
  void insitu_star_formation(struct galaxy_t* gal, int snapshot);

#ifdef __cplusplus
}
#endif

#endif
