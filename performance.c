/* cpropep.c  -  Calculation of Complex Chemical Equilibrium           */
/* Copyright (C) 2000                                                  */
/* Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                      */
/* Contribution of                                                     */
/* Mark Pinese <ida.pinese@bushnet.qld.edu.au>                         */

/* This program is free software; you can redistribute it and/or modify*/
/* it under the terms of the GNU General Public License as published by*/
/* the Free Software Foundation; either version 2 of the License, or   */
/* (at your option) any later version.                                 */

/* This program is distributed in the hope that it will be useful,     */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of      */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */ 
/* GNU General Public License for more details.                        */

/* You should have received a copy of the GNU General Public License   */
/* along with this program; if not, write to the Free Software         */
/* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.           */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <time.h>

#include "type.h"
#include "derivative.h"
#include "print.h"
#include "equilibrium.h"

#define TEMP_ITERATION_MAX  8

double velocity_of_flow(equilibrium_t *e, double exit_temperature);
double compute_temperature(equilibrium_t *e, double pressure,
                           double p_entropy);

double g = 9.80665;

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

double mixture_specific_heat_0(equilibrium_t *e, double temp)
{
  int i;
  double cp = 0.0;
  /* for gases */
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    cp += e->p->coef[GAS][i]*specific_heat_0(e->p->species[GAS][i], temp);
  }
  /* for condensed */
  for (i = 0; i < e->p->n[CONDENSED]; i++)
  {
    cp += e->p->coef[CONDENSED][i]*
      specific_heat_0(e->p->species[CONDENSED][i], temp);
  }
  return cp;
}

    
/* The temperature could be found by entropy conservation with a
   specified pressure */
double compute_temperature(equilibrium_t *e, double pressure, double p_entropy)
{
  int i = 0;

  double delta_lnt;
  //double p_entropy;   /* Enthalpy of the mixture (chamber) */
  double temperature;

  
  /* The first approximation is the chamber temperature */
  temperature = e->T; 
  //p_entropy = product_entropy(e);
  // e->entropy;//product_entropy(e);
  
  do
  {
    delta_lnt = (p_entropy  -
                 product_entropy_exit(e, pressure, temperature))
      /mixture_specific_heat_0(e, temperature);

    temperature = exp (log(temperature) + delta_lnt);

    
    //if (P == SP)
    //equilibrium(e, SP);
    
    
    i++;
    
  } while (fabs(delta_lnt) >= 0.5e-4 && i < TEMP_ITERATION_MAX);

  return temperature;
}

double velocity_of_flow(equilibrium_t *e, double exit_temperature)
{

  return sqrt(2000*(product_enthalpy(e)*R*e->T -
                 product_enthalpy_exit(e, exit_temperature)*
                 R*exit_temperature)/propellant_mass(e));

}

int frozen_performance(equilibrium_t *e, double exit_pressure)
{
  double sound_velocity;
  double flow_velocity;
  double throat_temperature;
  double exit_temperature;
  double pc_pt;
  double cp;
  double cp_cv; /* Ratio of specific heat, it is equal to the
                   isentropic exponent for frozen problem */
  double chamber_entropy;

  clock_t timer;
  timer = clock();

  
  printf("\nFrozen performance characteristics.\n");
  printf("-----------------------------------\n");
  
  /* Print the propellant composition */
  print_propellant_composition(e);

  /* find the equilibrium composition in the chamber */
  if (equilibrium(e, HP))
  {
    printf("No equilibrium, performance evaluation aborted.\n");
    return ERROR;
  }

  /* Cp of the combustion point assuming frozen */ 
  cp = mixture_specific_heat_0(e, e->T)*R;
  cp_cv = cp / (cp - e->n*R);
  

  printf("\nChamber results.\n");
  printf("-----------------\n");
  if (e->verbose > 0)
  {
    print_product_composition(e);  
  }
  print_product_properties(e);
  printf("Cp                   : % .2f\n", cp);
  printf("Cp/Cv                : % .2f\n", cp_cv);
    
  chamber_entropy = product_entropy(e);
  
  /* Find the exit temperature */
  exit_temperature = compute_temperature(e, exit_pressure, chamber_entropy);
  
  /* We must check if the exit temperature is more than 50 K lower
     than any transition temperature of condensed species */

  
  printf("\nExit conditions\n");
  printf("---------------\n");
  printf("Exit temperature  : % .2f K\n",  exit_temperature);
  printf("Exit pressure     : % .2f atm\n", exit_pressure);
  printf("Velocity of flow  : % .2f m/s\n",
         velocity_of_flow(e, exit_temperature));
  printf("Specific impulse  : % .2f s\n",
         velocity_of_flow(e, exit_temperature)/g);


  /* Computing throat condition */

  /* Approximation of the throat pressure */
  pc_pt = pow(cp_cv/2 + 0.5, cp_cv/(cp_cv - 1) );
  do
  {
    throat_temperature = compute_temperature(e, e->P/pc_pt, chamber_entropy);
    sound_velocity = sqrt( 1000*e->n*R*throat_temperature*
                           cp_cv/propellant_mass(e));
    flow_velocity = velocity_of_flow(e, throat_temperature);
    pc_pt = pc_pt / ( 1 + ((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                           /(1000*(cp_cv + 1)*e->n*R*throat_temperature
                             /propellant_mass(e))));
    
  } while (fabs( (pow(flow_velocity, 2) - pow(sound_velocity, 2))
                 /pow(flow_velocity, 2)) > 0.4e-4);

  printf("\nThroat conditions\n");
  printf("------------------\n");
  printf("Pc/Pt             : % .2f \n", pc_pt);
  printf("Throat pressure   : % .2f atm\n", e->P/pc_pt);
  printf("Throat temperature: % .2f K\n", throat_temperature);
  printf("Sound velocity    : % .2f m/s\n", sound_velocity);

  printf("\nComputation time for frozen: %f s\n",
         (double)(clock() - timer)/CLOCKS_PER_SEC);
  printf("---------------------------------------------\n");
  
  return SUCCESS;
  
}


int equilibrium_performance(equilibrium_t *e, double exit_pressure)
{

  double sound_velocity;
  double flow_velocity;
  double throat_temperature;
  double pc_pt;
  double chamber_entropy;
  
  deriv_t d;

  clock_t timer;
  
  /* Allocate a new equlibrium structure to hold information of the
     equilibrium at different point we consider */
  equilibrium_t *ne;
  ne = (equilibrium_t *)malloc(sizeof(equilibrium_t));

  timer = clock();

  printf("\nEquilibrium performance characteristics.\n");
  printf("-----------------------------------------.\n");
  
  /* Print the propellant composition */
  print_propellant_composition(e);
  
  /* find the equilibrium composition in the chamber */
  if (equilibrium(e, HP))
  {
    printf("No equilibrium, performance evaluation aborted.\n");
    return ERROR;
  }

  derivative(e, &d);
  
  printf("\nChamber results.\n");
  printf("-----------------\n");
  if (e->verbose > 0)
  {
    print_product_composition(e);
  }
  print_product_properties(e);
  print_derivative_results(d);
  
  *ne = *e;  
  chamber_entropy = product_entropy(e);

  
  /* Computing throat condition */
  /* Approximation of the throat pressure */
  pc_pt = pow(d.isex/2 + 0.5, d.isex/(d.isex - 1) );
  ne->entropy = chamber_entropy;
  do
  {
    ne->P = e->P/pc_pt;

    if (equilibrium(ne, SP))
    {
      printf("No equilibrium, performance evaluation aborted.\n");
      break;
    }
    if (e->verbose > 0)
    {
      print_product_composition(e);
    }
    
    derivative(e, &d);
    throat_temperature = ne->T;
    sound_velocity = sqrt (1000*ne->n*R*throat_temperature*
                           d.isex/propellant_mass(ne));
    
    flow_velocity = sqrt (2000*(product_enthalpy(e)*R*e->T -
                                product_enthalpy(ne)*R*ne->T)/
                          propellant_mass(ne));

    pc_pt = pc_pt / ( 1 + ((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                           /(1000*(d.isex + 1)*ne->n*R*throat_temperature
                             /propellant_mass(ne))));
    
  } while (fabs( (pow(flow_velocity, 2) - pow(sound_velocity, 2))
                 /pow(flow_velocity, 2)) > 0.4e-4);

  derivative(e, &d);
  
  printf("\nThroat conditions\n");
  printf("-----------------\n");
  printf("Pc/Pt                : % .2f\n", pc_pt);
  printf("Sound velocity       : % .2f m/s\n", sound_velocity);

  print_product_properties(ne);
  print_derivative_results(d);
  
  ne->entropy = chamber_entropy;
  ne->P = exit_pressure;

  //printf("Allo\n");
  if (equilibrium(ne, SP))
  {
    printf("No equilibrium, performance evaluation aborted.\n");
    return ERROR;
  }
  //printf("Allo\n");
  
  derivative(e, &d);
  
  printf("\nExit conditions\n");
  printf("-----------------\n");
 
  flow_velocity = sqrt(2000*(product_enthalpy(e)*R*e->T -
                             product_enthalpy(ne)*R*ne->T)/
                       propellant_mass(ne));
  
  printf("Velocity of flow     : % .2f m/s\n", flow_velocity);
  printf("Specific impulse     : % .2f s \n", flow_velocity/g);
  
  print_product_properties(ne);
  print_derivative_results(d);

  printf("\nComputation time for equilibrium: %f s\n",
         (double)(clock() - timer)/CLOCKS_PER_SEC);
  
  return 0;
}



