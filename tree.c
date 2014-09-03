/*
 * Copyright 2009 by Mathias Helsen en Kelly Beernaert
 * Last update: Thu Nov 26 2009
 */
#include "tree.h"

int printGnuplotCube(sim_ptr simulation, treeNode_ptr hdr)
{
  if(hdr != NULL){
   double x1[3] = {hdr->center[0] - hdr->width/2., hdr->center[1] - hdr->width/2., hdr->center[2] - hdr->width/2.};
   double x2[3] = {hdr->center[0] + hdr->width/2., hdr->center[1] - hdr->width/2., hdr->center[2] - hdr->width/2.};
   double x3[3] = {hdr->center[0] + hdr->width/2., hdr->center[1] - hdr->width/2., hdr->center[2] + hdr->width/2.};
   double x4[3] = {hdr->center[0] - hdr->width/2., hdr->center[1] - hdr->width/2., hdr->center[2] + hdr->width/2.};
   double x5[3] = {hdr->center[0] - hdr->width/2., hdr->center[1] + hdr->width/2., hdr->center[2] - hdr->width/2.};
   double x6[3] = {hdr->center[0] + hdr->width/2., hdr->center[1] + hdr->width/2., hdr->center[2] - hdr->width/2.};
   double x7[3] = {hdr->center[0] + hdr->width/2., hdr->center[1] + hdr->width/2., hdr->center[2] + hdr->width/2.};
   double x8[3] = {hdr->center[0] - hdr->width/2., hdr->center[1] + hdr->width/2., hdr->center[2] + hdr->width/2.};

  // print node bounding box data
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x4[0],x4[1],x4[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x1[0],x1[1],x1[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x2[0],x2[1],x2[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x3[0],x3[1],x3[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x4[0],x4[1],x4[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x8[0],x8[1],x8[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x7[0],x7[1],x7[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x6[0],x6[1],x6[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x5[0],x5[1],x5[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x8[0],x8[1],x8[2]);

	fprintf(simulation->nodefile, "\n\n");
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x1[0],x1[1],x1[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x5[0],x5[1],x5[2]);

	fprintf(simulation->nodefile, "\n\n");
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x2[0],x2[1],x2[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x6[0],x6[1],x6[2]);

	fprintf(simulation->nodefile, "\n\n");
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x3[0],x3[1],x3[2]);
	fprintf(simulation->nodefile, "%f\t%f\t%f\n", x7[0],x7[1],x7[2]);

	fprintf(simulation->nodefile, "\n\n");
  	
	int i;
	for(i=0;i<8;i++){
	    printGnuplotCube(simulation, hdr->children[i]);
	}
	return(0);
	}else{
		return(1);
	}

}

inline int
getOctant(treeNode_ptr parent, particle_ptr newparticle, sim_ptr simulation)
{
	//These will carry the boolean values of the octant (positive or negative).
	int x_oct, y_oct, z_oct, octant;
	//Set the correct flags according to the octant numbering scheme.
	//E.g. if the particle is larger than the center in the x direction
	//this boolean flags will be set to 0 and we will be able to only address 
	//the first 4.
	x_oct = (parent->center[0] > newparticle->pos[0]) ? 1 : 0;
	y_oct = (parent->center[1] > newparticle->pos[1]) ? 1 : 0;
	z_oct = (parent->center[2] > newparticle->pos[2]) ? 1 : 0;
	//Calculate the octant index.
	octant = z_oct + y_oct * 2 + x_oct * 4;

	return(octant);
}

int
deleteNode(treeNode_ptr parent, int octant, sim_ptr simulation)
{
	//If the to be deleted node exists...
	if(parent->children[octant] != NULL)
	{
		//delete its children using this method...
		int i;
		for(i=0;i<8;i++)
		{
			deleteNode(parent->children[octant], i, simulation);
		}
		//...and now delete the node itself.
		free(parent->children[octant]);
		parent->children[octant] = NULL;
		//Report succes
		return(0);
	}
	else
	{
		//If the to be deleted node is empty, something
		//might have gone wrong. Report this.
		return(1);
	}
}

int
addNode(treeNode_ptr hdr, particle_ptr newparticle, sim_ptr simulation)
{
	if(hdr->mass == 0)
	{
		//If the current node is empty, just fill in the correct parameters
		updateNode(hdr, newparticle, simulation);
	}
	else
	{
		if(hdr->noChildren == 1)
		{
			/*
			 * The current node is a leaf, therfore we will need to create a new
			 * subnode in the correct octant of the current node and give this 
			 * new subnode the parameters of the current node
			 */

			//Calculate the octant of the new node relative to this node
			int octant_subnode = getOctantCenterOfMass(hdr, simulation);
			//Create the node
			createNode(hdr, octant_subnode, simulation);
			//Shorthand notation for this new node
			treeNode_ptr newSubnode = hdr->children[octant_subnode];
			//Give this subnode the correct mass and center of mass location
			newSubnode->mass = hdr->mass;
			newSubnode->center_of_mass[0] = hdr->center_of_mass[0];
			newSubnode->center_of_mass[1] = hdr->center_of_mass[1];
			newSubnode->center_of_mass[2] = hdr->center_of_mass[2];
			//Ofcourse now this node has children
			hdr->noChildren = 0;
			//Try inserting the particle in this node again
			addNode(hdr, newparticle, simulation);
		}
		else
		{
			/*
			 * The current node is not a leaf, nor is it empty. Therefore
			 * we will update this node to contain the new particle and 
			 * try inserting the new particle into the correct octant
			 */
			//First update this node to contain the new particle
			updateNode(hdr, newparticle, simulation);
			//Calculate in which octant we want to insert the particle
			int octant = getOctant(hdr, newparticle, simulation);
			if(hdr->children[octant] == NULL)
			{
				/*
				 * The childnode in this octant does not exist,
				 * therefore we will create it and insert the 
				 * particle
				 */

				//Create the node
				createNode(hdr, octant, simulation);
				//Give it the parameters of the particle
				updateNode(hdr->children[octant], newparticle, simulation);
			}
			else
			{
				/*
				 * Now we have a non-null childnode, so we will
				 * try (using this method) to insert the particle
				 * in this node
				 */
				addNode(hdr->children[octant], newparticle, simulation);
			}
		}
	}
	return(0);
}

inline int
updateNode(treeNode_ptr parent, particle_ptr newparticle, sim_ptr simulation)
{
	//Shift the center of mass
	parent->center_of_mass[0] = parent->center_of_mass[0]*parent->mass
					+ newparticle->pos[0]*newparticle->mass;
	parent->center_of_mass[1] = parent->center_of_mass[1]*parent->mass
					+ newparticle->pos[1]*newparticle->mass;
	parent->center_of_mass[2] = parent->center_of_mass[2]*parent->mass
					+ newparticle->pos[2]*newparticle->mass;

	//Increase the total mass
	parent->mass += newparticle->mass;

	//Renormalize the center of mass
	parent->center_of_mass[0] /= parent->mass;
	parent->center_of_mass[1] /= parent->mass;
	parent->center_of_mass[2] /= parent->mass;

	return(0);
}

inline int
getOctantCenterOfMass(treeNode_ptr parent, sim_ptr simulation)
{
	//These will carry the boolean values of the octant (positive or negative).
	int x_oct, y_oct, z_oct, octant;
	//Set the correct flags according to the octant numbering scheme.
	//E.g. if the particle is larger than the center in the x direction
	//this boolean flags will be set to 0 and we will be able to only address 
	//the first 4.
	x_oct = (parent->center[0] > parent->center_of_mass[0]) ? 1 : 0;
	y_oct = (parent->center[1] > parent->center_of_mass[1]) ? 1 : 0;
	z_oct = (parent->center[2] > parent->center_of_mass[2]) ? 1 : 0;
	//Calculate the octant index.
	octant = z_oct + y_oct * 2 + x_oct * 4;

	return(octant);
	
}

int
createNode(treeNode_ptr parent, int octant, sim_ptr simulation)
{
  //check if the parent already has this subnode
  if(parent->children[octant] == NULL)
    {
      //we create the node
      treeNode_ptr node;
      node = malloc(sizeof(treeNode_type));
      parent->children[octant] = node;
      
      //calculate the coordinates of this node
      //These will carry the boolean values of the octant (positive or negative).
      int x_oct, y_oct, z_oct;
      //x_oct = (octant>=4) ? 1 : 0;
//      y_oct = ((octant - 4*x_oct)>=2) ? 1 : 0;
//      z_oct = octant - 4*x_oct - 2*y_oct;

      x_oct = (octant & 4)/4;
      y_oct = (octant & 2)/2;
      z_oct = (octant & 1);
      
      node->center[0] = parent->center[0] - (2.*(double)x_oct - 1)* parent->width * 0.25;
      node->center[1] = parent->center[1] - (2.*(double)y_oct - 1)* parent->width * 0.25;
      node->center[2] = parent->center[2] - (2.*(double)z_oct - 1)* parent->width * 0.25;
      
      //Set the width of the new node.
      node->width = parent->width * 0.5;
      
      //So far this node is childless.
      node->noChildren = 1;

      int i;
      for(i=0;i<8;i++)
      {
	node->children[i] = NULL;
      }

      //Initialise the mass of the node to zero;
      node->mass = 0;      

      //report succes.
      return(0);
    }
  else
    {
      //if this subnode already exits something may have gone wrong...
      return(1);
    }
}

int
makeTree(treeNode_ptr hdr, particlelist_ptr list, sim_ptr simulation)
{
	double width[3];
	width[0] = (list->xmax - list->xmin) * 1.005;
	width[1] = (list->ymax - list->ymin) * 1.005;
	width[2] = (list->zmax - list->zmin) * 1.005;
	double maxwidth;
	maxwidth = width[0];
	maxwidth = (maxwidth > width[1]) ? maxwidth : width[1];
	maxwidth = (maxwidth > width[2]) ? maxwidth : width[2];
	//The header's width is about the difference between the extrema of the coords
	hdr->width = maxwidth;
	//It has no mass 
	hdr->mass = 0;
	//No children either
	hdr->noChildren = 1;
	//It's centered around the width (i.e. symmetric)
	hdr->center[0] = 0.5*(list->xmax + list->xmin);
	hdr->center[1] = 0.5*(list->ymax + list->ymin);
	hdr->center[2] = 0.5*(list->zmax + list->zmin);
	//Center of mass is not yet applicable
	hdr->center_of_mass[0] = 0;
	hdr->center_of_mass[1] = 0;
	hdr->center_of_mass[2] = 0;

	int i;
	for(i=0; i<8; i++)
	{
		hdr->children[i] = NULL;
	}
	//For all the particles in the list, add to the tree
	//Parameter error checks if the addNode function was succesfull.
	int error=0;
	i=0;
	while(i<list->amount && error == 0)
	{
		error = addNode(hdr, list->deeltjeslijst[i], simulation);
		i++;
	}
	if(error == 0){
		return(0);
	}else{
		return(1);
	}
}

int
printLeaves(treeNode_ptr hdr, FILE * fileout, sim_ptr simulation)
{
	if(hdr == NULL){
		return(0);
	}else{
		if(hdr->noChildren == 1){
			fprintf(fileout, "%e\t%e\t%e\t%e\n",
					hdr->mass, hdr->center_of_mass[0], 
					hdr->center_of_mass[1], hdr->center_of_mass[2]);
		}
		else
		{
			int i;
			for(i=0;i<8;i++){
				printLeaves(hdr->children[i], fileout, simulation);
			}
		}
		return(0);
	}
}

int
printNodes(treeNode_ptr hdr, FILE * fileout, sim_ptr simulation)
{
	if(hdr != NULL)
	{
		fprintf(fileout, "%e\t%e\t%e\t%e\n",
			hdr->mass, hdr->center_of_mass[0], 
			hdr->center_of_mass[1], hdr->center_of_mass[2]);

		int i;
		for(i=0; i<8; i++)
		{
			printNodes(hdr->children[i], fileout, simulation);
		}
	}
	return(0);
}
