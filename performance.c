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

#include <stdio.h>
#include <math.h>

#include "equilibrium.h"

#define TEMP_ITERATION_MAX  8

double velocity_of_flow(equilibrium_t *e, double exit_temperature);
double compute_temperature(equilibrium_t *e, double pressure);

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

double product_specific_heat(equilibrium_t *e, double temp)
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
  
  /* find the equilibrium composition in the chamber */
  equilibrium(e, HP);
  
  /* Find the exit temperature */

 
  exit_temperature = compute_temperature(e, exit_pressure);
  
  /* We must check if the exit temperature is more than 50 K lower
     than any transition temperature of condensed species */

  printf("\nFrozen performance characteristics.\n");
  printf("-----------------------------------\n");
  printf("Exit temperature: %f K\n", exit_temperature);
  printf("Exit pressure   : %f atm\n", exit_pressure);
  printf("Velocity of flow: %f m/s\n", velocity_of_flow(e, exit_temperature));
  printf("Specific impulse: %f s\n", velocity_of_flow(e, exit_temperature)/g);


  /* Computing throat condition */

  /* Cp of the combustion point */
  cp = product_specific_heat(e, e->T)*R;
  cp_cv = cp / (cp - e->n*R);
 
  printf("Cp/Cv : %f\n", cp_cv);

  /* Approximation of the throat pressure */
  pc_pt = pow(cp_cv/2 + 0.5, cp_cv/(cp_cv - 1) );

  do
  {
    //printf("Pc/Pt: %f\n", pc_pt);

    throat_temperature = compute_temperature(e, e->P/pc_pt);

    sound_velocity = sqrt( 1000*e->n*R*throat_temperature*
                           cp_cv/propellant_mass(e));
    
    flow_velocity = velocity_of_flow(e, throat_temperature);

    //printf("Flow velocity: %f\n", flow_velocity);
    //printf("Sound velocity: %f\n", sound_velocity);
    
    pc_pt = pc_pt / ( 1 + ((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                           /(1000*(cp_cv + 1)*e->n*R*throat_temperature
                             /propellant_mass(e))));
    
  } while (fabs( (pow(flow_velocity, 2) - pow(sound_velocity, 2))
                 /pow(flow_velocity, 2)) > 0.4e-4);
  
  printf("Pc/Pt: %f\n", pc_pt);
  printf("Throat temperature: %f\n", throat_temperature);
  printf("Sound velocity: %f\n", sound_velocity);
  
  return 0;
  
}

double compute_temperature(equilibrium_t *e, double pressure)
{
  int i = 0;

  double delta_lnt;
  double p_entropy;   /* Enthalpy of the mixture (chamber) */
  double temperature;

  
  /* The first approximation is the chamber temperature */
  temperature = e->T; 
  p_entropy = product_entropy(e);
  
  do
  {
    delta_lnt = (p_entropy  -
                 product_entropy_exit(e, pressure, temperature))
      /product_specific_heat(e, temperature);

    temperature = exp (log(temperature) + delta_lnt);
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





