/*
 * Copyright 2009 by Mathias Helsen and Kelly Beernaert
 * Last update Sat Nov 7 2009
 */
#ifndef _PARTICLELIST_H
#define _PARTICLELIST_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "particle.h"
#include "sim.h"

//We maken een struct (verzameling van variabelen in een object) en
//noemen dit particle (type = struct particle)
struct particlelist
{
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
  int amount;
  particle_ptr * deeltjeslijst;
};

//We maken ook een afkorting van een pointer naar zulk een object
typedef struct particlelist * particlelist_ptr;
//We maken ook een afkorting van zulk een object (i.p.v. struct particle nu particle_type
typedef struct particlelist particlelist_type;

/*
 *Makes a list of particles using random number generators
 *
 *This function makes a list of N particles with coördinates inside a cube with boundary coordinates -maxC and maxC
 *
 *@param list:         The list of particles to be made
 *@param N:            The number of particles in the list
 *@param maxC:         The maximum coordinate of the particles
 */
int MakeList(particlelist_ptr list, int N, double maxC, sim_ptr simulation);

/*
 *Deletes the given list
 *
 *@param list:          The list of particles.
 */

int DeleteList(particlelist_ptr list);

/*
 * Prints information about the particles in the list
 *
 * This function prints all the components of position and velocity
 * and mass of the particles in the specified list.
 *
 * @param list:		The list of particles.
 * @param fileout:	The file to which the data is to be printed, NULL prints to stdout.
 * @return		0 for succes, 1 for failure
 */
int printList(particlelist_ptr list, FILE * fileout, sim_ptr simulation);

/*
 * Reads a file containing particles and makes a particlelist from it
 *
 * This function does exactly the opposite that function printList does.
 * Files written in the format of the latter function are read and
 * a particlelist is created from it.
 * WARNING: Method does not check for array bounds!
 */
int parseList(particlelist_ptr list, FILE * filein, sim_ptr simulation);

#endif
