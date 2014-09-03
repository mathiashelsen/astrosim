/*
 * Copyright 2009 by Mathias Helsen en Kelly Beernaert
 * Last update: Thu Nov 26 2009
 */

#ifndef _ALGO_H
#define _ALGO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "particle.h"
#include "tree.h"
#include "sim.h"

/*
 *Function to calculate the gravitational force of all the other particles in the tree on ONE specified particle
 *
 *@param hdr:         The header of the tree.
 *@param particle:    The particle on which the force acts
 *@param force[]:     a 3D array to hold the resulting gravitational force of the calculation
 *@param trigger[]:   a 3D array with conditions for which we need to dig deeper into the tree
 */
int getForce(treeNode_ptr hdr, particle_ptr particle, double force[], sim_ptr simulation);

/*
 *The initialisation of the LeapFrog integrator
 */
inline int initLeapFrog(treeNode_ptr hdr, particlelist_ptr particles, sim_ptr simulation);

/*
 *Executes a single step in the LeapFrog integration scheme
 */
inline int leapFrogStep(treeNode_ptr hdr, particlelist_ptr particles, sim_ptr simulation);

/*
 *Solves the equation of motion for all the particles using a LeapFrog integration scheme
 */
int leapFrog(particlelist_ptr particles, sim_ptr simulation);

/*
 * Function calculates the potential energy FOR ONE PARTICLE using a tree method. Just like
 * the getForce function, this does not calculate it in N^2, but in in NlogN time.
 * However is uses an approximation.
 *
 * @param hdr	The pointer to the header of the tree
 * @param particles	The list of particles
 * @param trigger	The trigger for seeing if we need to go down into the node
 * @return	The total potential energy
 */
double getPotentialEnergy(treeNode_ptr hdr, particle_ptr particle, sim_ptr simulation);

/*
 * This function uses getPotentialEnergy to calculate the total potential energy
 */
double getTotalPotentialEnergy(treeNode_ptr hdr, particlelist_ptr particles, sim_ptr simulation);

/*
 * This function calculates the total kinetic energy by iterating over all the particles
 * in the particlelist.
 *
 * @param particles	The list of particles.
 * @return	The total kinetic energy
 */
inline double getTotalKineticEnergy(particlelist_ptr particles);

#endif
