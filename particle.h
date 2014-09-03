/*
 * Copyright 2009 by Mathias Helsen en Kelly Beernaert
 * Last update: Sat Nov 7 2009
 */

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <stdlib.h>
#include <stdio.h>

//We maken een struct (verzameling van variabelen in een object) en
//noemen dit particle (type = struct particle)
struct particle
{
  double mass;
  double pos[3];
  double vel[3];
};

//We maken ook een afkorting van een pointer naar zulk een object
typedef struct particle * particle_ptr;
//We maken ook een afkorting van zulk een object (i.p.v. struct particle nu particle_type
typedef struct particle particle_type;

#endif
