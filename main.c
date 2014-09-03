/*
 * Copyright 2009 by Mathias Helsen en Kelly Beernaert
 * Last update: Sat Nov 7 2009
*/


#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "particlelist.h"
#include "tree.h"
#include "sim.h"
#include "algo.h"

void test1 (double stepsize);
void test2 (void);
void test3 (int numparts);


int
main(int argc, char **argv)
{	
	if(argc==1){
		printf("./astroproject steps stepsPerPrint\n");
	}else{
		test3(atoi(argv[1]));
	}
	return 0;
}

void test3(int numparts){
		sim_ptr simulation;
		simulation = malloc(sizeof(sim_type));
		simulation->steps = 10;
		simulation->stepsPerPrint = 1;
		simulation->timestep = 0.01;
		simulation->logfile = fopen("logfile.log", "w");
		simulation->nodefile = fopen("nodefile.dat", "w");
		simulation->particlelistfile = fopen("particles.dat", "w");
		simulation->energyfile = fopen("energy.dat", "w");
		simulation->G = 1.0;
		simulation->trigger = 0.6;
		simulation->softening = 0.002;
		
		particlelist_ptr list = malloc(sizeof(particlelist_type));
		MakeList(list, numparts, 1.0, simulation);
		
		//treeNode_ptr hdr = malloc(sizeof(treeNode_type));
		//makeTree(hdr, list, simulation);
		//printGnuplotCube(simulation, hdr);
		//int i;
		//for(i=0;i<8;i++){
		//	deleteNode(hdr, i, simulation);
		//}
		
		leapFrog(list, simulation);
		
		//free(hdr);
		DeleteList(list);
		fclose(simulation->logfile);
		fclose(simulation->nodefile);
		fclose(simulation->particlelistfile);
		fclose(simulation->energyfile);
		free(simulation);
}
		
void test1(double stepsize){
		sim_ptr simulation;
		simulation = malloc(sizeof(sim_type));
		simulation->steps = 1000;
		simulation->stepsPerPrint = 1;
		simulation->timestep = stepsize;
		simulation->logfile = fopen("logfile.log", "w");
		simulation->nodefile = fopen("nodefile.dat", "w");
		simulation->particlelistfile = fopen("particles.dat", "w");
		simulation->energyfile = fopen("energy.dat", "w");
		simulation->G = 1.0;
		simulation->trigger = 0.6;
		simulation->softening = 0.0;

		particlelist_ptr particles = malloc(sizeof(particlelist_type));
		particles->deeltjeslijst = malloc(2*sizeof(particle_ptr));
		particles->amount = 2;
		particles->deeltjeslijst[0] = malloc(sizeof(particle_type));
		particles->deeltjeslijst[1] = malloc(sizeof(particle_type));

		particles->deeltjeslijst[0]->mass = 1.0;
		particles->deeltjeslijst[0]->pos[0] = 0.0;
		particles->deeltjeslijst[0]->pos[1] = 0.0;
		particles->deeltjeslijst[0]->pos[2] = 0.0;
		particles->deeltjeslijst[0]->vel[0] = -0.5;
		particles->deeltjeslijst[0]->vel[1] = 0.0;
		particles->deeltjeslijst[0]->vel[2] = 0.0;
	
		particles->deeltjeslijst[1]->mass = 1.0;
		particles->deeltjeslijst[1]->pos[0] = 0.0;
		particles->deeltjeslijst[1]->pos[1] = 1.0;
		particles->deeltjeslijst[1]->pos[2] = 0.0;
		particles->deeltjeslijst[1]->vel[0] = 0.5;
		particles->deeltjeslijst[1]->vel[1] = 0.0;
		particles->deeltjeslijst[1]->vel[2] = 0.0;

		particles->xmin = particles->xmax = particles->zmax = particles->zmin = particles->ymin = 0.0;
		particles->ymax = 1.0;

		treeNode_ptr header = malloc(sizeof(treeNode_type));
		makeTree(header, particles, simulation);
		//Test the force calculation
		double force[3] = {0.0, 0.0, 0.0};
		getForce(header, particles->deeltjeslijst[0], force, simulation);
		printf("Force1 = {%f, %f, %f}\n", force[0], force[1], force[2]);
		force[0] = force[1] = force[2] = 0.0;
		getForce(header, particles->deeltjeslijst[1], force, simulation);
		printf("Force2 = {%f, %f, %f}\n", force[0], force[1], force[2]);
		//Test the energy calculation
		double potential_energy = getTotalPotentialEnergy(header, particles, simulation);
		double kinetic_energy = getTotalKineticEnergy(particles);
		printf("Kinetic energy: %f\nPotential energy: %f\nTotal energy: %f\n", kinetic_energy,
				potential_energy, potential_energy + kinetic_energy);
		//Test the leapfrog integrator
		leapFrog(particles, simulation);

		//Start clearing the memory
		int i;
		for(i=0;i<8;i++){
			deleteNode(header, i, simulation);
		}
		free(header);
		DeleteList(particles);
		fclose(simulation->logfile);
		fclose(simulation->nodefile);
		fclose(simulation->particlelistfile);
		fclose(simulation->energyfile);
		free(simulation);
}

void test2(void){
		sim_ptr simulation;
		simulation = malloc(sizeof(sim_type));
		simulation->steps = 10;
		simulation->stepsPerPrint = 1;
		simulation->timestep = 0.0001;
		simulation->logfile = fopen("logfile.log", "w");
		simulation->nodefile = fopen("nodefile.dat", "w");
		simulation->particlelistfile = fopen("particles.dat", "w");
		simulation->energyfile = fopen("energy.dat", "w");
		simulation->G = 1.0;
		simulation->trigger = 0.6;
		simulation->softening = 0.0;

		particlelist_ptr particles = malloc(sizeof(particlelist_type));
		particles->deeltjeslijst = malloc(3*sizeof(particle_ptr));
		particles->amount = 3;
		particles->deeltjeslijst[0] = malloc(sizeof(particle_type));
		particles->deeltjeslijst[1] = malloc(sizeof(particle_type));
		particles->deeltjeslijst[2] = malloc(sizeof(particle_type));

		particles->deeltjeslijst[0]->mass = 1.0;
		particles->deeltjeslijst[0]->pos[0] = 0.0;
		particles->deeltjeslijst[0]->pos[1] = 0.0;
		particles->deeltjeslijst[0]->pos[2] = 0.0;
		particles->deeltjeslijst[0]->vel[0] = 0.0;
		particles->deeltjeslijst[0]->vel[1] = 0.0;
		particles->deeltjeslijst[0]->vel[2] = 0.0;
	
		particles->deeltjeslijst[1]->mass = 1e6;
		particles->deeltjeslijst[1]->pos[0] = 1000.0;
		particles->deeltjeslijst[1]->pos[1] = 0.0;
		particles->deeltjeslijst[1]->pos[2] = 1.0;
		particles->deeltjeslijst[1]->vel[0] = 0.0;
		particles->deeltjeslijst[1]->vel[1] = 0.0;
		particles->deeltjeslijst[1]->vel[2] = 0.0;

		particles->deeltjeslijst[2]->mass = 1e6;
		particles->deeltjeslijst[2]->pos[0] = 1000.0;
		particles->deeltjeslijst[2]->pos[1] = 0.0;
		particles->deeltjeslijst[2]->pos[2] = -1.0;
		particles->deeltjeslijst[2]->vel[0] = 0.0;
		particles->deeltjeslijst[2]->vel[1] = 0.0;
		particles->deeltjeslijst[2]->vel[2] = 0.0;
		
		particles->xmin = 0.0;
		particles->xmax = 1000.0;
		particles->zmax = 1.0;
		particles->zmin = -1.0;
		particles->ymin = 0.0;
		particles->ymax = 0.0;



		treeNode_ptr header = malloc(sizeof(treeNode_type));
		makeTree(header, particles, simulation);
		//Print a nice tree
		printGnuplotCube(simulation, header);
		//Test the force calculation
		double force[3] = {0.0, 0.0, 0.0};
		getForce(header, particles->deeltjeslijst[0], force, simulation);
		printf("Acc1 = {%f, %f, %f}\n", force[0], force[1], force[2]);
		force[0] = force[1] = force[2] = 0.0;
		getForce(header, particles->deeltjeslijst[1], force, simulation);
		printf("Acc2 = {%f, %f, %f}\n", force[0], force[1], force[2]);
		force[0] = force[1] = force[2] = 0.0;
		getForce(header, particles->deeltjeslijst[2], force, simulation);
		printf("Acc3 = {%f, %f, %f}\n", force[0], force[1], force[2]);
		//Test the energy calculation
		double potential_energy = getTotalPotentialEnergy(header, particles, simulation);
		double kinetic_energy = getTotalKineticEnergy(particles);
		printf("Kinetic energy: %f\nPotential energy: %f\nTotal energy: %f\n", kinetic_energy,
				potential_energy, potential_energy + kinetic_energy);
		//Test the leapfrog integrator
		leapFrog(particles, simulation);

		//Start clearing the memory
		int i;
		for(i=0;i<8;i++){
			deleteNode(header, i, simulation);
		}
		free(header);
		DeleteList(particles);
		fclose(simulation->logfile);
		fclose(simulation->nodefile);
		fclose(simulation->particlelistfile);
		fclose(simulation->energyfile);
		free(simulation);
}
