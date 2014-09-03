/*
 * Copyright 2009 by Mathias Helsen en Kelly Beernaert
 * Last update: Thu Nov 26 2009
 */
#ifndef _TREE_H
#define _TREE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "particle.h"
#include "particlelist.h"
#include "sim.h"

/* This struct contains all the components of one of the nodes of the tree. */
struct treeNode
{
	//The children of the node
	struct treeNode * children[8];
	//The center of the node
	double center[3];
	//The center of mass of the node
	double center_of_mass[3];
	//The width of the node
	double width;
	//The total mass included in the node and its subnodes
	double mass;

	// If node contains no children 1, otherwise 0
	int noChildren;
};

/* This typedef makes a type treeNode_ptr from a pointer to a struct treeNode */
typedef struct treeNode * treeNode_ptr;
/* This typedef makes a type treeNode_type from a struct treeNode */
typedef struct treeNode treeNode_type;

/*
 *Function graphically draws the tree.
 *
 *@param simulation:    The simulation for which the tree is printed
 *@param hdr:           The tree one wishes to print
 *@return               0 for succes, 1 for failure
 */
int printGnuplotCube(sim_ptr simulation, treeNode_ptr hdr); 

/*
 * Function creates one node.
 *
 * This function creates one node at the specified address (i.e. octant
 * of the parent node). Initializes all values at zero or true.
 *
 * @param parent:	The parent of the new node.
 * @param octant:	The octant in which the node is to be placed.
 * @return		0 for succes, 1 for failure.
 */
int createNode(treeNode_ptr parent, int octant, sim_ptr simulation);

/*
 * Function adds one node to the tree.
 *
 * This is a recursive function which finds the correct location of a particle
 * in the tree and then insert a new node (fresh made) containing all data 
 * about a particle.
 * As the code proceeds through the tree, it updates the included mass,
 * center of mass, ... to include the new particle
 *
 * @param hdr:		A treeNode_ptr to the header of the tree.
 * @param newparticle:	The particle that one wants to insert.
 * @return		0 for succes, 1 for failure.
 */
int addNode(treeNode_ptr hdr, particle_ptr newparticle, sim_ptr simulation);

/*
 * Function deletes a node and all of its subnodes in a tree.
 *
 * Function deletes a certain specified node in a certain octant
 * of the parent node and deletes all the subnodes in the process.
 *
 * @param parent:	The node that contains the to be deleted subnode.
 * @param octant:	The octant in which the to be deleted subnode is located
 * @return		0 for succes, 1 for failure.
 */
int deleteNode(treeNode_ptr parent, int octant, sim_ptr simulation);

/*
 * Function calculates the index number of the correct octant in which to place
 * the new particle.
 *
 * The function calculates the octant index using a scheme where octant 0 corresponds
 * with x,y,z positive, 1 with x,y positive and z negative, 2 with x and z positive 
 * and y negative, 3 with x positive and y and z negative, ..., 7 with x,y,z all
 * negative. This labeling scheme translates to a certain binary like octant labeling.
 *
 * @param parent:	The treeNode_ptr to the parent of the octants from which to choose.
 * @param newparticle:	The particle which is to be inserted.
 * @return		The required octant	
 */
inline int getOctant(treeNode_ptr parent, particle_ptr newparticle, sim_ptr simulation);

/*
 * Function calculates the octant index of the center of mass of a certain node.
 *
 * This function calculates the octant index of the center of mass of a certain
 * node. This function is required when making a subnode from a node (i.e. 
 * make a single particle node into a multiple particle node).
 * 
 * @param parent:	The node of which this needs to calculate the octant
 * @return		The required octant
 */
inline int getOctantCenterOfMass(treeNode_ptr parent, sim_ptr simulation);

/*
 * Function updates parameters of a node when a particle is added.
 *
 * This function calculates the new parameters (e.g. mass) of a node when
 * a new particle is added to this node or one if its subnodes.
 *
 * @param parent:	The node which is to be updated.
 * @param newparticle:	The particle which is inserted.
 * @return		0 for succes, 1 for failure
 */
inline int updateNode(treeNode_ptr parent, particle_ptr newparticle, sim_ptr simulation);

/*
 * Function creates a tree from a particlelist
 *
 * This function takes a particlelist and creates a tree from it.
 *
 * @param hdr:		The header of the tree.
 * @param list:		The list of the particles.
 * @return		0 for succes, 1 for failure.
 */
int makeTree(treeNode_ptr hdr, particlelist_ptr list, sim_ptr simulation);

/*
 *Function prints the particles and their properties
 *
 *@param hdr:           The header of the tree.
 *@param fileout:       The name of the file to print to.
 */
int printLeaves(treeNode_ptr hdr, FILE * fileout, sim_ptr simulation);

/*
 *Function prints the nodes and their properties
 *
 *@param hdr:           The header of the tree.
 *@param fileout:       The name of the file to print to.
 */
int printNodes(treeNode_ptr hdr, FILE * fileout, sim_ptr simulation);

#endif
