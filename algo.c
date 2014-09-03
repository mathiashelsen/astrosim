/*
 * Copyright 2009 by Mathias Helsen en Kelly Beernaert
 * Last update: Thu Nov 26 2009
 */
#include "algo.h"

int getForce(treeNode_ptr hdr, particle_ptr particle, double force[], sim_ptr simulation)
{
  double dif[3], dist;
  
  int i;
  for(i=0; i<8; i++)
    {
      //check if the node exists
      if(hdr->children[i] != NULL)
	{
	  dif[0] = particle->pos[0] - hdr->children[i]->center_of_mass[0];
	  dif[1] = particle->pos[1] - hdr->children[i]->center_of_mass[1];
	  dif[2] = particle->pos[2] - hdr->children[i]->center_of_mass[2];
	  dist = sqrt(dif[0]*dif[0] + dif[1]*dif[1] + dif[2]*dif[2]);
		  
	  if(dist > 0.00001) //we do not want to calculate the force of the particle on the particle itself
	    { 
		if((hdr->children[i]->width/dist > simulation->trigger) && hdr->children[i]->noChildren == 0)
		{
		  //we need to dig deeper in the tree
	//	  printf("Go deeper into the tree.\n");
		  getForce(hdr->children[i], particle, force, simulation);
		}
	      else
		{
	//	if(hdr->children[i]->noChildren == 1){
	//		printf("This node has no children.\n");
	//	}else{
	//		printf("We accept the approximation.\n");
	//	}
		double tempvar;
		tempvar = dist*dist + simulation->softening*simulation->softening;
		tempvar *= sqrt(tempvar);
		  force[0]+= -(simulation->G)*hdr->children[i]->mass*dif[0]/tempvar;
		  force[1]+= -(simulation->G)*hdr->children[i]->mass*dif[1]/tempvar;
		  force[2]+= -(simulation->G)*hdr->children[i]->mass*dif[2]/tempvar;
		}
	    }
	}
    }
  return(0);
}

double
getPotentialEnergy(treeNode_ptr hdr, particle_ptr particle, sim_ptr simulation)
{
	int i;
	double energy=0;
	for(i=0;i<8;i++){
		if(hdr->children[i] != NULL){
			double separation[3], distance;

			separation[0] = particle->pos[0] - hdr->children[i]->center_of_mass[0];
			separation[1] = particle->pos[1] - hdr->children[i]->center_of_mass[1];
			separation[2] = particle->pos[2] - hdr->children[i]->center_of_mass[2];
			distance = sqrt(separation[0]*separation[0] + separation[1]*separation[1]
				+ separation[2]*separation[2]);

			if(distance > 0.00001){
				if(hdr->children[i]->width/distance > simulation->trigger && hdr->children[i]->noChildren==0){
					energy += getPotentialEnergy(hdr->children[i], particle, simulation);
				}else{
					double tempvar;
					tempvar = sqrt(distance*distance + simulation->softening*simulation->softening);
					energy -= (simulation->G) * hdr->children[i]->mass * particle->mass / tempvar;
				}
			}
		}
	}
	return(energy);
}

double
getTotalPotentialEnergy(treeNode_ptr hdr, particlelist_ptr particles, sim_ptr simulation)
{
	double energy=0;
	int i, amount = particles->amount;
	for(i=0;i<amount;i++){
		energy += getPotentialEnergy(hdr, particles->deeltjeslijst[i], simulation);
	}
	return(energy/2.0);
}

inline double
getTotalKineticEnergy(particlelist_ptr particles)
{
	double energy=0;
	int i, amount = particles->amount;
	particle_ptr particle;
	for(i=0;i<amount;i++){
		particle = particles->deeltjeslijst[i];
		energy += particle->mass * 0.5 * (particle->vel[0]*particle->vel[0]
			+ particle->vel[1]*particle->vel[1]
			+ particle->vel[2]*particle->vel[2]);
	}
	return(energy);
}
	
inline int 
initLeapFrog(treeNode_ptr hdr, particlelist_ptr particles, sim_ptr simulation)
{
  particle_ptr particle;
  int i, amount = particles->amount;
  double force[3], timestep = simulation->timestep;
  for(i=0;i<amount;i++)
    {
      particle = particles->deeltjeslijst[i];
      force[0] = force[1] = force[2] = 0.;
      getForce(hdr, particle,force, simulation);
      particle->pos[0] = particle->pos[0] 
	+ particle->vel[0]*timestep/2.
	+ force[0]*timestep*timestep/8.; 
      particle->pos[1] = particle->pos[1] 
	+ particle->vel[1]*timestep/2.
	+ force[1]*timestep*timestep/8.; 
      particle->pos[2] = particle->pos[2] 
	+ particle->vel[2]*timestep/2.
	+ force[2]*timestep*timestep/8.; 
    }
  return(0);
}

inline int leapFrogStep(treeNode_ptr hdr, particlelist_ptr particles, sim_ptr simulation)
{
  particle_ptr particle;
  int i, amount = particles->amount;
  double force[3], timestep = simulation->timestep;
  makeTree(hdr, particles, simulation);
	particles->xmax = particles->xmin = particles->deeltjeslijst[0]->pos[0];
	particles->ymax = particles->ymin = particles->deeltjeslijst[0]->pos[1];
	particles->zmax = particles->zmin = particles->deeltjeslijst[0]->pos[2];
  for(i=0;i<amount;i++)
    {
      particle = particles->deeltjeslijst[i];
      force[0] = force[1] = force[2] = 0.;
      getForce(hdr, particle,force, simulation);
      particle->vel[0] = particle->vel[0] + force[0]*timestep;
      particle->vel[1] = particle->vel[1] + force[1]*timestep;
      particle->vel[2] = particle->vel[2] + force[2]*timestep;
      particle->pos[0] = particle->pos[0] + particle->vel[0]*timestep; 
      particle->pos[1] = particle->pos[1] + particle->vel[1]*timestep; 
      particle->pos[2] = particle->pos[2] + particle->vel[2]*timestep;

      particles->xmin = (particle->pos[0] < particles->xmin) ? particle->pos[0] : particles->xmin;
      particles->xmax = (particle->pos[0] > particles->xmax) ? particle->pos[0] : particles->xmax;
      particles->ymin = (particle->pos[1] < particles->ymin) ? particle->pos[1] : particles->ymin;
      particles->ymax = (particle->pos[1] > particles->ymax) ? particle->pos[1] : particles->ymax;
      particles->zmin = (particle->pos[2] < particles->zmin) ? particle->pos[2] : particles->zmin;
      particles->zmax = (particle->pos[2] > particles->zmax) ? particle->pos[2] : particles->zmax;
    }
  for(i=0;i<8;i++)
    {
      deleteNode(hdr, i, simulation);
    }
  
  return(0);
}

int leapFrog(particlelist_ptr particles, sim_ptr simulation)
{
  //Variables for holding the energies
  double kinetic_energy, potential_energy, total_energy, delta, total_energy_init;
  //Allocate a header for the tree
  treeNode_ptr hdr = malloc(sizeof(treeNode_type));
  //Make a tree
  makeTree(hdr, particles, simulation);
  //Calculate the different energies
  kinetic_energy = getTotalKineticEnergy(particles);
  potential_energy = getTotalPotentialEnergy(hdr, particles, simulation);
  total_energy_init = kinetic_energy + potential_energy;

  //Initialize the leapfrog integrator
  initLeapFrog(hdr, particles, simulation);

  //Tear down the tree
  int i,j,k;
  for(i=0; i<8; i++){
	deleteNode(hdr, i, simulation);
  }

  for(i=0; i<(simulation->steps / simulation->stepsPerPrint); i++)
    {
      for(j=0;j<simulation->stepsPerPrint;j++){
      	leapFrogStep(hdr, particles, simulation);
      }
	  //Print out the position of the particles and their velocity
     // char filename[40];
     // FILE * fileout;
     // sprintf(filename, "particles_%d.dat", i+1);
     // fileout = fopen(filename, "w");
	  printList(particles, simulation->particlelistfile, simulation);
//	  fclose(fileout);
	  fprintf(simulation->particlelistfile, "\n\n");
	  //Create a new tree
  	  makeTree(hdr, particles, simulation);
	  //Calculate the energies
	  kinetic_energy = getTotalKineticEnergy(particles);
	  potential_energy = getTotalPotentialEnergy(hdr, particles, simulation);
	  total_energy = kinetic_energy + potential_energy;
	  //Calculate the relative difference in energy (compared to the start)
	  delta = fabs((total_energy_init - total_energy) / total_energy_init);
	  fprintf(simulation->energyfile, "%e\t%e\t%e\t%e\n", kinetic_energy, potential_energy,
			  total_energy, delta);
	  fflush(simulation->energyfile);
	  //Tear down the tree
	 for(k=0;k<8;k++){
		 deleteNode(hdr, k, simulation);
	 }
    }
  free(hdr);
  return(0);
}
