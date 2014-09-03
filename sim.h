/*
 * Copyright 2009 by Mathias Helsen and Kelly Beernaert
 * Last update Thu Nov 26 2009
 */
#ifndef _SIM_H
#define _SIM_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct sim
{
  //steps:		the number of integration step
  //stepsPerPrint:	the number of integration steps performed before writting the current particlelist to a file
  int steps, stepsPerPrint;
  //the timestep of the integration
  double timestep;
  //the file to which to log
  FILE * logfile;
  //the file in which all the nodes are printed
  FILE * nodefile;

  FILE * particlelistfile;

  FILE * energyfile;

  double G, trigger, softening;
  
  //several time keeping variabels.
  struct tm *local;
  time_t t;
};

typedef struct sim sim_type;
typedef struct sim * sim_ptr;

#endif
