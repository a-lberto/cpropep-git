/* performance.c  -  Compute performance caracteristic of a motor
                     considering equilibrium                      */
/* $Id: performance.c,v 1.8 2000/06/14 00:27:50 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese  <pinese@cyberwizards.com.au>                        */
/*                                                                     */
/* Licensed under the GPLv2                                            */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <time.h>

#include "performance.h"
#include "derivative.h"
#include "print.h"
#include "equilibrium.h"

#include "compat.h"
#include "return.h"

#define TEMP_ITERATION_MAX  8

double velocity_of_flow(equilibrium_t *e, double exit_temperature);
double compute_temperature(equilibrium_t *e, double pressure,
                           double p_entropy);

double g = 9.80665; /* m/s/s */

extern FILE * errorfile;
extern FILE * outputfile;

/* Entropy of the product at the exit pressure and temperature */
double product_entropy_exit(equilibrium_t *e, double pressure, double temp)
{
  double ent;
  double t = e->T;
  double p = e->P;
  e->T = temp;
  e->P = pressure;
  ent = product_entropy(e); /* entropy at the new pressure */
  e->T = t;
  e->P = p;
  return ent;
}

/* Enthalpy of the product at the exit temperature */
double product_enthalpy_exit(equilibrium_t *e, double temp)
{
  double enth;
  double t = e->T;
  e->T = temp;
  enth = product_enthalpy(e);
  e->T = t;
  return enth;
}

/* The specific heat of the mixture for frozen performance */
double mixture_specific_heat_0(equilibrium_t *e, double temp)
{
  int i;
  double cp = 0.0;
  /* for gases */
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    cp += e->p.coef[i][GAS]*specific_heat_0(e->p.species[i][GAS], temp);
  }
  /* for condensed */
  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {
    cp += e->p.coef[i][CONDENSED]*
      specific_heat_0(e->p.species[i][CONDENSED], temp);
  }
  return cp;
}

    
/* The temperature could be found by entropy conservation with a
   specified pressure */
double compute_temperature(equilibrium_t *e, double pressure, double p_entropy)
{
  int i = 0;

  double delta_lnt;
  double temperature;

  /* The first approximation is the chamber temperature */
  temperature = e->T;
  
  do
  {
    delta_lnt = (p_entropy  -
                 product_entropy_exit(e, pressure, temperature))
      /mixture_specific_heat_0(e, temperature);

    temperature = exp (log(temperature) + delta_lnt);
        
    i++;

  } while (fabs(delta_lnt) >= 0.5e-4 && i < TEMP_ITERATION_MAX);

  return temperature;
}

double velocity_of_flow(equilibrium_t *e, double exit_temperature)
{
  return sqrt(2000*(product_enthalpy(e)*R*e->T -
                    product_enthalpy_exit(e, exit_temperature)*
                    R*exit_temperature));//  /propellant_mass(e));
}

int frozen_performance(equilibrium_t *e, performance_t *p,
                       double exit_pressure)
{
  double sound_velocity;
  double flow_velocity;
  double pc_pt;
  double chamber_entropy;
  

  /* find the equilibrium composition in the chamber */
  if (!(e->isequil)) /* if the equilibrium have not already been compute */
  {
    if (equilibrium(e, HP))
    {
      fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
      return ERROR;
    }
  }

  p->frozen.chamber.temperature = e->T;
  p->frozen.chamber.pressure    = e->P;
  p->frozen.chamber.velocity    = 0.0;
  
  /* Cp of the combustion point assuming frozen */
  p->frozen.cp    = mixture_specific_heat_0(e, e->T)*R;
  /* Cv = Cp - nR  (for frozen) */
  p->frozen.cp_cv = p->frozen.cp / (p->frozen.cp - e->n*R);

  if (e->verbose > 0)
  {
    print_product_composition(e);  
  }

  p->frozen.molar_mass = product_molar_mass(e);
  chamber_entropy      = product_entropy(e);
  
  /* Find the exit temperature */
  p->frozen.exit.temperature = compute_temperature(e, exit_pressure,
                                                   chamber_entropy);
  
  /* We must check if the exit temperature is more than 50 K lower
     than any transition temperature of condensed species.
     In this case the results are not good and must be reject. */

  p->frozen.exit.pressure    = exit_pressure;
  p->frozen.exit.velocity    = velocity_of_flow(e, p->frozen.exit.temperature);
  p->frozen.specific_impulse = p->frozen.exit.velocity/g;

  /* aera per unit mass flow rate */
  p->frozen.exit.aera_dotm   = 1000 * R * p->frozen.exit.temperature * e->n /
    (p->frozen.exit.pressure * p->frozen.exit.velocity);

  
  
  /* Computing throat condition */
  pc_pt = pow (p->frozen.cp_cv/2 + 0.5,
               p->frozen.cp_cv/(p->frozen.cp_cv - 1) );
  do
  {
    p->frozen.throat.temperature = compute_temperature(e, e->P/pc_pt,
                                                       chamber_entropy);
    
    sound_velocity = sqrt( 1000 * e->n * R * p->frozen.throat.temperature *
                           p->frozen.cp_cv);
    
    flow_velocity = velocity_of_flow(e, p->frozen.throat.temperature);
    
    pc_pt = pc_pt / ( 1 + ((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                           /(1000*(p->frozen.cp_cv + 1)*
                             e->n*R*p->frozen.throat.temperature)));
    
  } while (fabs( (pow(flow_velocity, 2) - pow(sound_velocity, 2))
                 /pow(flow_velocity, 2)) > 0.4e-4);


  p->frozen.throat.pressure    = e->P/pc_pt;
  p->frozen.throat.velocity    = sound_velocity;

  /* aera per unit mass flow rate */
  p->frozen.throat.aera_dotm   = 1000 * R * p->frozen.throat.temperature
    * e->n /
    (p->frozen.throat.pressure * p->frozen.throat.velocity);

  
  p->frozen_ok = true;
  
  return SUCCESS;
}


int equilibrium_performance(equilibrium_t *e, equilibrium_t *ne,
                            performance_t *p, double exit_pressure)
{

  double sound_velocity = 0.0;
  double flow_velocity;
  double pc_pt;
  double chamber_entropy;
  
  deriv_t *d;
  
  /* Allocate a new equlibrium structure to hold information of the
     equilibrium at different point we consider */

  //equilibrium_t *ne;
  //ne = (equilibrium_t *)malloc(sizeof(equilibrium_t));
  //initialize_equilibrium(ne);

  /* find the equilibrium composition in the chamber */
  if (!(e->isequil))/* if the equilibrium have not already been compute */
  {
    if (equilibrium(e, HP))
    {
      fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
      return ERROR;
    }
  }
  p->equilibrium.chamber.pressure = e->P;
  p->equilibrium.chamber.temperature = e->T;
  p->equilibrium.chamber.molar_mass = product_molar_mass(e);
  
  /* first consider the chamber state */
  d = &(p->equilibrium.chamber.deriv);
  
  /* Compute derivative */
  derivative(e, d);

  if (e->verbose > 0)
  {
    print_product_composition(e);
  }

  /* Begin by first aproximate the new equilibrium to be
     the same as the chamber equilibrium */

  copy_equilibrium(ne, e);
  
  chamber_entropy = product_entropy(e);
  
  /* Computing throat condition */
  /* Approximation of the throat pressure */
  pc_pt = pow(d->isex/2 + 0.5, d->isex/(d->isex - 1) );
  ne->entropy = chamber_entropy;
  do
  {
    ne->P = e->P/pc_pt;

    /* We must compute the new equilibrium each time */
    if (equilibrium(ne, SP))
    {
      fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
      return ERROR;
    }

    if (e->verbose > 0)
    {
      print_product_composition(e);
    }

    /* now consider the throat state */
    d = &(p->equilibrium.throat.deriv);
    
    /* Compute the new derivative properties */
    derivative(ne, d);
    p->equilibrium.throat.temperature = ne->T;

    sound_velocity = sqrt (1000*ne->n*R*p->equilibrium.throat.temperature*
                           d->isex);///propellant_mass(ne));
    
    flow_velocity = sqrt (2000*(product_enthalpy(e)*R*e->T -
                                product_enthalpy(ne)*R*ne->T)
                          );///propellant_mass(ne));

    pc_pt = pc_pt / ( 1 + ((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                           /(1000*(d->isex + 1)*ne->n*R*
                             p->equilibrium.throat.temperature
                             )));///propellant_mass(ne))));
    
  } while (fabs( (pow(flow_velocity, 2) - pow(sound_velocity, 2))
                 /pow(flow_velocity, 2)) > 0.4e-4);

  derivative(ne, d);

  p->equilibrium.throat.pressure   = e->P/pc_pt;
  p->equilibrium.throat.velocity   = sound_velocity;
  p->equilibrium.throat.molar_mass = product_molar_mass(ne);

  p->equilibrium.throat.aera_dotm   = 1000 * R *
    p->equilibrium.throat.temperature
    * e->n /
    (p->equilibrium.throat.pressure * p->equilibrium.throat.velocity);
  
  ne->entropy = chamber_entropy;
  ne->P = exit_pressure;

  /* Find the exit equilibrium */
  if (equilibrium(ne, SP))
  {
    fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
    return ERROR;
  }

  /* now consider the exit state */
  d = &(p->equilibrium.exit.deriv); 
  
  /* Compute new derivative properties */
  derivative(ne, d);

  p->equilibrium.exit.molar_mass = product_molar_mass(ne);
 
  flow_velocity = sqrt(2000*(product_enthalpy(e)*R*e->T -
                             product_enthalpy(ne)*R*ne->T)
                       );///propellant_mass(ne));

  p->equilibrium.exit.temperature = ne->T;
  p->equilibrium.exit.pressure    = ne->P;
  p->equilibrium.exit.velocity    = flow_velocity;
  p->equilibrium.specific_impulse = flow_velocity/g;

  p->equilibrium.exit.aera_dotm   = 1000 * R *
    p->equilibrium.exit.temperature
    * e->n /
    (p->equilibrium.exit.pressure * p->equilibrium.exit.velocity);
  
  p->equilibrium_ok = true;

  /* we have to deallocate the new equilibrium */
  //dealloc_equilibrium (ne);
  //free (ne);
  
  return SUCCESS;
}





