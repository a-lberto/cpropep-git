/* derivative.c  -  Fill the mattrix to compute thermochemical derivative
                    relative to logarithm of pressure and temperature */
/* $Id: derivative.c,v 1.7 2000/06/14 00:27:50 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <pinese@cyberwizards.com.au>                         */
/*                                                                     */
/* Licensed under the GPLv2                                            */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "derivative.h"
#include "libnum.h"
#include "equilibrium.h"
#include "performance.h"

#include "compat.h"
#include "return.h"

extern FILE * errorfile;
extern FILE * outputfile;

/* Compute the specific_heat of the mixture using thermodynamics
   derivative with respect to logarithm of temperature */
double mixture_specific_heat(equilibrium_t *e, double *sol)
{
  int i, j;
  double cp, tmp;
  
  cp = 0.0;
  /* Compute Cp/R */
  for (i = 0; i < e->n_element; i++)
  {
    tmp = 0.0;
    for (j = 0; j < e->p.n[GAS]; j++)
      tmp += product_element_coef (e->element[i], e->p.species[j][GAS])*
        e->p.coef[j][GAS]*enthalpy_0(e->p.species[j][GAS], e->T);
    
    cp += tmp*sol[i];
    
  }
  
  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {
    cp += enthalpy_0(e->p.species[i][CONDENSED], e->T)
      *sol[i + e->n_element];
  }
  
  tmp = 0.0;
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    tmp += e->p.coef[i][GAS] * enthalpy_0(e->p.species[i][GAS], e->T);
  }
  cp += tmp*sol[e->n_element + e->p.n[CONDENSED]];
  
  cp += mixture_specific_heat_0(e, e->T);
  
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    cp += e->p.coef[i][GAS]*pow(enthalpy_0(e->p.species[i][GAS], e->T),2);
  }
  //for (i = 0; i < e->p.n[CONDENSED]; i++)
  //{
  //  cp += e->p.coef[i][CONDENSED]*
  //    pow(enthalpy_0(e->p.species[i][CONDENSED], e->T),2);
  // }
  return cp;
}

int derivative(equilibrium_t *e, deriv_t *d)
{
//  int i;
  int size;
#ifdef TRUE_ARRAY
  double *matrix;
#else
  double ** matrix;   
#endif
  double  * sol;

  /* the size of the coefficient matrix */
  size = e->n_element + e->p.n[CONDENSED] + 1;

#ifdef TRUE_ARRAY
  matrix = (double *) malloc (size * (size + 1) * sizeof(double));
#else
  /* allocate the memory for the matrix */
  matrix = (double **) malloc (sizeof(double *) * size);
  for (i = 0; i < size; i++)
    matrix[i] = (double *) malloc (sizeof(double) * (size+1));
#endif
  
  
  /* allocate the memory for the solution vector */
  sol = (double *) calloc (size, sizeof(double));

  fill_temperature_derivative_matrix(matrix, e);
  
  if ( lu(matrix, sol, size) == -1 )
  {
    fprintf(outputfile, "The matrix is singular.\n");
  }
  else
  {
    if (e->verbose > 1)
    {
      fprintf(outputfile, "Temperature derivative results.\n");
      print_vec(sol, size);
    }
    
    d->cp = mixture_specific_heat(e, sol)*R;
    d->del_lnV_lnT = 1 + sol[e->n_element + e->p.n[CONDENSED]];
      
  }

  fill_pressure_derivative_matrix(matrix, e);

  if ( lu(matrix, sol, size) == -1 )
  {
    fprintf(outputfile, "The matrix is singular.\n");
  }
  else
  {
    if (e->verbose > 1)
    {
      fprintf(outputfile, "Pressure derivative results.\n");
      print_vec(sol, size);
    }
    d->del_lnV_lnP = sol[e->n_element + e->p.n[CONDENSED]] - 1;
    
  }

  d->cv    = d->cp + e->n*R * pow(d->del_lnV_lnT, 2)/d->del_lnV_lnP;
  d->cp_cv = d->cp/d->cv;
  d->isex  = -d->cp_cv / d->del_lnV_lnP;
  d->vson  = sqrt(1000 * e->n * R * e->T * d->isex);
  
//  if (e->verbose > 1)
//  {
//    fprintf(outputfile, "\n");
//    fprintf(outputfile, "del ln(V)/del ln(T)  = % f\n", d->del_lnV_lnT);
//    fprintf(outputfile, "del ln(V)/del ln(P)  = % f\n", d->del_lnV_lnP);
//    fprintf(outputfile, "Cp                   = % f\n", d->cp);
//    fprintf(outputfile, "Cv                   = % f\n", d->cv);
//    fprintf(outputfile, "Cp/Cv                = % f\n", d->cp_cv);
//    fprintf(outputfile, "Isentropic exponent  = % f\n", d->isex);
//    fprintf(outputfile, "RT/V                 = % f\n", e->P/e->n);
//  }
  
#ifndef TRUE_ARRAY
  for (i = 0; i < size; i++)
    free(matrix[i]);
#endif
  
  free(matrix);
  free(sol);
  return 0;
}


/* Fill the matrix with the coefficient for evaluating derivatives with
   respect to logarithm of temperature at constant pressure */
#ifdef TRUE_ARRAY
int fill_temperature_derivative_matrix(double *matrix, equilibrium_t *e)
#else
int fill_temperature_derivative_matrix(double **matrix, equilibrium_t *e)
#endif
{
  
  int j, k, size;
  double tmp;

  /* the size of the coefficient matrix */
  size = e->n_element + e->p.n[CONDENSED] + 1;
  
  /* fill the common part */
  fill_matrix(matrix, e, TP);
  
  /* del ln(n)/ del ln(T) */
#ifdef TRUE_ARRAY
  *(matrix + e->n_element + e->p.n[CONDENSED] + size *
    (e->n_element + e->p.n[CONDENSED])) = 0.0;
#else
  matrix[e->n_element + e->p.n[CONDENSED]][e->n_element + e->p.n[CONDENSED]]
    =  0.0;
#endif
  
  /* right side */
  for (j = 0; j < e->n_element; j++)
  {
    tmp = 0.0;
    for (k = 0; k < e->p.n[GAS]; k++)
      tmp -= product_element_coef( e->element[j], e->p.species[k][GAS]) * 
        e->p.coef[k][GAS] * enthalpy_0(e->p.species[k][GAS], e->T);
#ifdef TRUE_ARRAY
    *(matrix + j + size * (e->n_element + e->p.n[CONDENSED] + 1)) = tmp;
#else
    matrix[j][ e->n_element + e->p.n[CONDENSED] + 1] = tmp;
#endif
  }

  for (j = 0; j < e->p.n[CONDENSED]; j++) /* row */
#ifdef TRUE_ARRAY
    *(matrix + j + e->n_element + size *
      (e->n_element + e->p.n[CONDENSED] + 1)) =
      -enthalpy_0( e->p.species[j][CONDENSED], e->T);
#else
    matrix[ j + e->n_element ][ e->n_element + e->p.n[CONDENSED] + 1] = 
      -enthalpy_0( e->p.species[j][CONDENSED], e->T);
#endif
  
  tmp = 0.0;
  for (k = 0; k < e->p.n[GAS]; k++)
    tmp -= e->p.coef[k][GAS]*enthalpy_0( e->p.species[k][GAS], e->T ); 

#ifdef TRUE_ARRAY
  *(matrix + e->n_element + e->p.n[CONDENSED] + size *
    (e->n_element + e->p.n[CONDENSED]+ 1)) = tmp;
#else
  matrix[e->n_element + e->p.n[CONDENSED]][e->n_element
                                           + e->p.n[CONDENSED]+ 1] = tmp;
#endif
  
  
  return 0;
}

/* Fill the matrix with the coefficient for evaluating derivatives with
   respect to logarithm of pressure at constant temperature */
#ifdef TRUE_ARRAY
int fill_pressure_derivative_matrix(double *matrix, equilibrium_t *e)
#else
int fill_pressure_derivative_matrix(double **matrix, equilibrium_t *e)
#endif
{
  
  int j, k, size;
  double tmp;

  /* the size of the coefficient matrix */
  size = e->n_element + e->p.n[CONDENSED] + 1;
  
  /* fill the common part */
  fill_matrix(matrix, e, TP);
  
  /* del ln(n)/ del ln(T) */
#ifdef TRUE_ARRAY
  *(matrix + e->n_element + e->p.n[CONDENSED] + size *
    (e->n_element + e->p.n[CONDENSED])) = 0.0;
#else
  matrix[e->n_element + e->p.n[CONDENSED]][e->n_element + e->p.n[CONDENSED]]
    =  0.0;
#endif
  
  
  /* right side */
  for (j = 0; j < e->n_element; j++)
  {
    tmp = 0.0;
    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += product_element_coef( e->element[j], e->p.species[k][GAS]) * 
        e->p.coef[k][GAS];
#ifdef TRUE_ARRAY
    *(matrix + j + size * (e->n_element + e->p.n[CONDENSED] + 1)) = tmp;
#else
    matrix[j][ e->n_element + e->p.n[CONDENSED] + 1] = tmp;
#endif
  }

  
  for (j = 0; j < e->p.n[CONDENSED]; j++) /* row */
#ifdef TRUE_ARRAY
    *(matrix + j + e->n_element + size *
      (e->n_element + e->p.n[CONDENSED] + 1)) = 0.0;
#else
    matrix[ j + e->n_element ][ e->n_element + e->p.n[CONDENSED] + 1] = 0;
#endif
  
  tmp = 0.0;
  for (k = 0; k < e->p.n[GAS]; k++)
    tmp += e->p.coef[k][GAS]; 

#ifdef TRUE_ARRAY
  *(matrix + e->n_element + e->p.n[CONDENSED] + size *
    (e->n_element + e->p.n[CONDENSED]+ 1)) = tmp;
#else
  matrix[e->n_element + e->p.n[CONDENSED]][e->n_element
                                          + e->p.n[CONDENSED]+ 1] = tmp;
#endif
  
  return 0;
}

