/*
 * Copyright 2009 by Mathias Helsen and Kelly Beernaert
 * Last update Sat Nov 7 2009
 */
#include "particlelist.h"
#include <gsl/gsl_rng.h>

#define LINESIZE 1024

int MakeList(particlelist_ptr list, int N, double maxC, sim_ptr simulation)
{
  list->deeltjeslijst = malloc(sizeof(particle_ptr)*N);
  if(list->deeltjeslijst == NULL)
  {
    return(1);
  }else{
    list->xmin = maxC;
    list->ymin = maxC;
    list->zmin = maxC;
    list->xmax = -maxC;
    list->ymax = -maxC;
    list->zmax = -maxC;
  
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    double x,y,z;
    int i;
    for(i=0; i<N; i++)
      {
        list->deeltjeslijst[i]=malloc(sizeof(particle_type));
	if(list->deeltjeslijst[i] == NULL){
		return(1);
	}else{
	        x = (gsl_rng_uniform_pos(r) - 0.5)*2*maxC;
        	y = (gsl_rng_uniform_pos(r) - 0.5)*2*maxC;
	        z = (gsl_rng_uniform_pos(r) - 0.5)*2*maxC;
        	list->xmin = (x < list->xmin) ? x : list->xmin;
	        list->xmax = (x > list->xmax) ? x : list->xmax;
	        list->ymin = (y < list->ymin) ? y : list->ymin;
        	list->ymax = (y > list->ymax) ? y : list->ymax;
	        list->zmin = (z < list->zmin) ? z : list->zmin;
	        list->zmax = (z > list->zmax) ? z : list->zmax;
       
        	list->deeltjeslijst[i]->pos[0]=x;
	        list->deeltjeslijst[i]->pos[1]=y;
        	list->deeltjeslijst[i]->pos[2]=z;
	        list->deeltjeslijst[i]->vel[0]= gsl_rng_uniform(r);
	        list->deeltjeslijst[i]->vel[1]= gsl_rng_uniform(r);
	        list->deeltjeslijst[i]->vel[2]= gsl_rng_uniform(r);
		list->deeltjeslijst[i]->mass = 1e-30;
	}
      }
    list->amount = N;
    
    gsl_rng_free(r); 
    return(0);
  }
}
  
int
printList(particlelist_ptr list, FILE * fileout, sim_ptr simulation)
{
	if(fileout == NULL) fileout = stdout;
	if(list == NULL)
	{
		//List does not exist, return error
		return(1);
	}
	else
	{
		//Print the header
		fprintf(fileout, "#id\tmass\t\tx coord\t\ty coord\t\tz coord\t\tx vel\t\ty vel\t\tz vel\n");
		int i;
		for(i=0;i<list->amount;i++)
		{
			//Print the parameters of all the particles
			fprintf(fileout, "%d\t%e\t%f\t%f\t%f\t%f\t%f\t%f\n", i, list->deeltjeslijst[i]->mass, 
					list->deeltjeslijst[i]->pos[0], list->deeltjeslijst[i]->pos[1], 
					list->deeltjeslijst[i]->pos[2], list->deeltjeslijst[i]->vel[0],
					list->deeltjeslijst[i]->vel[1], list->deeltjeslijst[i]->vel[2]);
		}
		//Operation succesfull
		return(0);
	}
}

//TODO: CHECK FOR EXTENTS OF PARTICLES
int
parseList(particlelist_ptr list, FILE * filein, sim_ptr simulation)
{
	double xmax=0.0, xmin=0.0, ymax=0.0, ymin=0.0, zmax=0.0, zmin=0.0;
	if(filein == NULL){
		return(1);
	}else{
		//We need word and line string
		char *word, *line;
		int i=0;
		//Allocate for a line (maximum chars per line is LINESIZE)
		line = malloc(LINESIZE * sizeof(char));
		//Read (and skip) the first line
		if(fgets(line, LINESIZE, filein) == NULL){
			return(1);
		}
		//Chop the file into lines
		while(fgets(line, LINESIZE, filein) != NULL){
			int j=0;
			//Chop the line with a TAB as separator
			word = strtok(line, "\t");
			while(word!=NULL){
				switch(j){
					//Firstly allocate a particle
					//Set this particle's mass
					case 0: list->deeltjeslijst[i] = malloc(sizeof(particle_type));
						list->deeltjeslijst[i]->mass = atof(word);
						break;

					//Set the particle's location
					case 1: list->deeltjeslijst[i]->pos[0] = atof(word);
						xmax = (list->deeltjeslijst[i]->pos[0] > xmax) ? list->deeltjeslijst[i]->pos[0] : xmax;
						xmin = (list->deeltjeslijst[i]->pos[0] < xmin) ? list->deeltjeslijst[i]->pos[0] : xmin;
						break;
	
					case 2: list->deeltjeslijst[i]->pos[1] = atof(word);
						ymax = (list->deeltjeslijst[i]->pos[1] > ymax) ? list->deeltjeslijst[i]->pos[1] : ymax;
						ymin = (list->deeltjeslijst[i]->pos[1] < ymin) ? list->deeltjeslijst[i]->pos[1] : ymin;
						break;
	
					case 3: list->deeltjeslijst[i]->pos[2] = atof(word);
						zmax = (list->deeltjeslijst[i]->pos[2] > zmax) ? list->deeltjeslijst[i]->pos[2] : zmax;
						zmin = (list->deeltjeslijst[i]->pos[2] < zmin) ? list->deeltjeslijst[i]->pos[2] : zmin;
						break;

					//And set the particle's velocity
					case 4: list->deeltjeslijst[i]->vel[0] = atof(word);
						break;

					case 5: list->deeltjeslijst[i]->vel[1] = atof(word);
						break;

					case 6: list->deeltjeslijst[i]->vel[2] = atof(word);
						break;
				}
				//Chop again
				word = strtok(NULL, "\t");
				j++;
			}
			i++;
		}
		//Ofcourse set the number of particles we have added
		list->amount = i;
		//We where succesfull
		return(0);

	}
}

int DeleteList(particlelist_ptr list)
{
  int i;
  for(i=0; i<list->amount; i++)
    {
      free(list->deeltjeslijst[i]);
    }
  free(list->deeltjeslijst);
  free(list);
  return(0);
}


