/* performance.c  -  Compute performance caracteristic of a motor
                     considering equilibrium                      */
/* $Id: performance.c,v 1.9 2000/06/20 02:15:12 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese  <pinese@cyberwizards.com.au>                        */
/*                                                                     */
/* Licensed under the GPLv2                                            */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "performance.h"
#include "derivative.h"
#include "print.h"
#include "equilibrium.h"

#include "const.h"
#include "compat.h"
#include "return.h"
#include "thermo.h"

#define TEMP_ITERATION_MAX  8

double velocity_of_flow(equilibrium_t *e, double exit_temperature);
double compute_temperature(equilibrium_t *e, double pressure,
                           double p_entropy);


/* Entropy of the product at the exit pressure and temperature */
double product_entropy_exit(equilibrium_t *e, double pressure, double temp)
{
  double ent;
  double t = e->properties.T;
  double p = e->properties.P;
  e->properties.T = temp;
  e->properties.P = pressure;
  ent = product_entropy(e); /* entropy at the new pressure */
  e->properties.T = t;
  e->properties.P = p;
  return ent;
}

/* Enthalpy of the product at the exit temperature */
double product_enthalpy_exit(equilibrium_t *e, double temp)
{
  double enth;
  double t = e->properties.T;
  e->properties.T = temp;
  enth = product_enthalpy(e);
  e->properties.T = t;
  return enth;
}
    
/* The temperature could be found by entropy conservation with a
   specified pressure */
double compute_temperature(equilibrium_t *e, double pressure, double p_entropy)
{
  int i = 0;

  double delta_lnt;
  double temperature;

  /* The first approximation is the chamber temperature */
  temperature = e->properties.T;
  
  do
  {
    delta_lnt = (p_entropy  - product_entropy_exit(e, pressure, temperature))
      / mixture_specific_heat_0(e, temperature);

    temperature = exp (log(temperature) + delta_lnt);
        
    i++;

  } while (fabs(delta_lnt) >= 0.5e-4 && i < TEMP_ITERATION_MAX);

  return temperature;
}

double velocity_of_flow(equilibrium_t *e, double exit_temperature)
{
  return sqrt(2000*(product_enthalpy(e)*R*e->properties.T -
                    product_enthalpy_exit(e, exit_temperature)*
                    R*exit_temperature));
}

int frozen_performance(equilibrium_t *e, double exit_pressure)
{
  double sound_velocity;
  double flow_velocity;
  double pc_pt;
  double cp_cv;
  double chamber_entropy;
  
  equilibrium_t *t  = e + 1;  /* throat equilibrium */
  equilibrium_t *ex = e + 2; /* exit equilibrium   */
  
  /* find the equilibrium composition in the chamber */
  if (!(e->product.isequil))
  {
    if (equilibrium(e, HP))
    {
      fprintf(outputfile,
              "No equilibrium, performance evaluation aborted.\n");
      return ERROR;
    }
  }

  e->properties.dV_T =  1.0;
  e->properties.dV_P = -1.0;
  e->properties.Cp   = mixture_specific_heat_0(e, e->properties.T) * R;
  e->properties.Cv   = e->properties.Cp - e->itn.n * R;
  e->properties.Isex = e->properties.Cp/e->properties.Cv;

  compute_thermo_properties(e);
  
  chamber_entropy  = product_entropy(e);
  
  /* begin computation of throat caracteristic */
  copy_equilibrium(t, e);

  t->properties.dV_T =  1.0;
  t->properties.dV_P = -1.0;
  
  cp_cv = e->properties.Cp/e->properties.Cv;
  /* first estimate of Pc/Pt */
  pc_pt = pow (cp_cv/2 + 0.5, cp_cv/(cp_cv - 1));

  do
  {
   
    t->properties.T = compute_temperature(t, e->properties.P/pc_pt,
                                          chamber_entropy);

    /* Cp of the combustion point assuming frozen */
    t->properties.Cp   = mixture_specific_heat_0(e, t->properties.T) * R;
    /* Cv = Cp - nR  (for frozen) */
    t->properties.Cv   = t->properties.Cp - t->itn.n * R;
    t->properties.Isex = t->properties.Cp/t->properties.Cv;

    compute_thermo_properties(t);
    
    sound_velocity = sqrt(1000 * e->itn.n * R * t->properties.T *
                          t->properties.Isex);
    
    flow_velocity = velocity_of_flow(e, t->properties.T);
    
    pc_pt = pc_pt / ( 1 + ((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                           /(1000*(t->properties.Isex + 1)*
                             t->itn.n * R *t->properties.T)));
    
  } while (fabs( (pow(flow_velocity, 2) - pow(sound_velocity, 2))
                 /pow(flow_velocity, 2)) > 0.4e-4);


  t->properties.P    = e->properties.P/pc_pt;
  t->performance.Isp = t->properties.Vson = sound_velocity;

  /* Now compute exit properties */ 
  copy_equilibrium(ex, e);
  
  ex->properties.T = compute_temperature(e, exit_pressure,
                                         chamber_entropy);
  /* We must check if the exit temperature is more than 50 K lower
     than any transition temperature of condensed species.
     In this case the results are not good and must be reject. */

  ex->properties.P       = exit_pressure;
  ex->performance.Isp    = velocity_of_flow(e, ex->properties.T);
  ex->performance.a_dotm = 1000 * R * ex->properties.T * ex->itn.n /
    (ex->properties.P * ex->performance.Isp);

  ex->properties.dV_T =  1.0;
  ex->properties.dV_P = -1.0;
  /* Cp of the combustion point assuming frozen */
  ex->properties.Cp   = mixture_specific_heat_0(e, ex->properties.T) * R;
  /* Cv = Cp - nR  (for frozen) */
  ex->properties.Cv   = ex->properties.Cp - ex->itn.n * R;
  ex->properties.Isex = ex->properties.Cp/ex->properties.Cv;

  compute_thermo_properties(ex);
  
  ex->properties.Vson = sqrt( 1000 * e->itn.n * R * ex->properties.T *
                              e->properties.Isex);
    
  t->performance.a_dotm = 1000 * R * t->properties.T
    * t->itn.n / (t->properties.P * t->performance.Isp);
  

  t->performance.ae_at = 1.0;
  t->performance.cstar = e->properties.P * t->performance.a_dotm;
  t->performance.cf    = t->performance.Isp /
    (e->properties.P * t->performance.a_dotm);
  t->performance.Ivac  = t->performance.Isp + t->properties.P
    * t->performance.a_dotm;

  ex->performance.ae_at =
    (ex->properties.T * t->properties.P * t->performance.Isp) /
    (t->properties.T * ex->properties.P * ex->performance.Isp);
  ex->performance.cstar = e->properties.P * t->performance.a_dotm;
  ex->performance.cf    = ex->performance.Isp /
    (e->properties.P * t->performance.a_dotm);
  ex->performance.Ivac  = ex->performance.Isp + ex->properties.P
    * ex->performance.a_dotm;

  return SUCCESS;
}


int equilibrium_performance(equilibrium_t *e, double exit_pressure)
{

  double sound_velocity = 0.0;
  double flow_velocity;
  double pc_pt;
  double chamber_entropy;
  
  equilibrium_t *t  = e + 1; // throat equilibrium
  equilibrium_t *ex = e + 2; // throat equilibrium
  
  /* find the equilibrium composition in the chamber */
  if (!(e->product.isequil))
/* if the equilibrium have not already been compute */
  {
    if (equilibrium(e, HP))
    {
      fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
      return ERROR;
    }
  }

  /* Begin by first aproximate the new equilibrium to be
     the same as the chamber equilibrium */

  copy_equilibrium(t, e);
  
  chamber_entropy = product_entropy(e);
  
  /* Computing throat condition */
  /* Approximation of the throat pressure */
  pc_pt = pow(t->properties.Isex/2 + 0.5,
              t->properties.Isex/(t->properties.Isex - 1) );

  t->entropy = chamber_entropy;

  do
  {
    t->properties.P = e->properties.P/pc_pt;

    /* We must compute the new equilibrium each time */
    if (equilibrium(t, SP))
    {
      fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
      return ERROR;
    }

    sound_velocity = sqrt (1000*t->itn.n*R*t->properties.T*
                           t->properties.Isex);
    
    flow_velocity = sqrt (2000*(product_enthalpy(e)*R*e->properties.T -
                                product_enthalpy(t)*R*t->properties.T));

    pc_pt = pc_pt / ( 1 + ((pow(flow_velocity, 2) - pow(sound_velocity, 2))
                           /(1000*(t->properties.Isex + 1)*t->itn.n*R*
                             t->properties.T)));
    
  } while (fabs( (pow(flow_velocity, 2) - pow(sound_velocity, 2))
                 /pow(flow_velocity, 2)) > 0.4e-4);


  t->properties.P = e->properties.P/pc_pt;
  t->properties.Vson = sound_velocity;
  t->performance.Isp = sound_velocity;


  t->performance.a_dotm = 1000 * R *
    t->properties.T * t->itn.n /
    (t->properties.P * t->performance.Isp);
  
  copy_equilibrium(ex, e);

  ex->entropy      = chamber_entropy;
  ex->properties.P = exit_pressure;

  /* Find the exit equilibrium */
  if (equilibrium(ex, SP))
  {
    fprintf(outputfile, "No equilibrium, performance evaluation aborted.\n");
    return ERROR;
  }
  
  flow_velocity = sqrt(2000*(product_enthalpy(e)*R*e->properties.T -
                             product_enthalpy(ex)*R*ex->properties.T));

  
  ex->performance.Isp = flow_velocity;

  ex->performance.a_dotm = 1000 * R *
    ex->properties.T * ex->itn.n /
    (ex->properties.P * ex->performance.Isp);
 
  t->performance.ae_at = 1.0;
  t->performance.cstar = e->properties.P * t->performance.a_dotm;
  t->performance.cf    = t->performance.Isp /
    (e->properties.P * t->performance.a_dotm);
  t->performance.Ivac  = t->performance.Isp + t->properties.P
    * t->performance.a_dotm;
  
  ex->performance.ae_at =
    (ex->properties.T * t->properties.P * t->performance.Isp) /
    (t->properties.T * ex->properties.P * ex->performance.Isp);
  ex->performance.cstar = e->properties.P * t->performance.a_dotm;
  ex->performance.cf    = ex->performance.Isp /
    (e->properties.P * t->performance.a_dotm);
  ex->performance.Ivac  = ex->performance.Isp + ex->properties.P
    * ex->performance.a_dotm;
  
  return SUCCESS;
}





