/* num.c  -  Library of numerical method                               */
/* Copyright (C) 2000                                                  */
/* Antoine Lefebvre <antoine.lefebvre@polymtl.ca                       */

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
#include <stdlib.h>
#include <math.h>

//#include <time.h>

#include "num.h"

FILE * errorfile;
FILE * outputfile;

int function(int neq, double time, double *y, double *dy, 
             void *data);

double f(double x);

double f1(double x);
double df1(double x);

double g1(double x);

int main(void)
{

  int i, n;
  double *ans;
  double *ic;  

  double *matrix;
  double *solution;

  int size = 3;

  printf("Testing the LU factorisation algotythm.\n");
  matrix = (double *) malloc (sizeof(double)*size*(size+1));
  solution = (double *) malloc (sizeof(double)*size);
  
  matrix[0] = 0;
  matrix[1] = 1;
  matrix[2] = 3;
  matrix[3] = 2;
  matrix[4] = 0;
  matrix[5] = 0;
  matrix[6] = 1;
  matrix[7] = 0;
  matrix[8] = 1;
  matrix[9] = -2/3;
  matrix[10] = 5/2;
  matrix[11] = 1;

  //NUM_matscale(matrix, size);
  NUM_print_matrix(matrix, size);

  NUM_lu(matrix, solution, size);

  NUM_print_vec(solution, size);

/*
  printf("Secante\n");
  printf("Solution: %f\n",  sec(f, -15, 0, 100, 0.0001));

  printf("Newton\n");
  printf("Solution: %f\n", newton(f1, df1, 1, 100, 0.0001));
  
  printf("Ptfix\n");
  printf("Solution: %f\n", ptfix(g1, 1, 100, 0.0001));

  //printf("epsilon: %.16e\n", epsilon());
*/  

  printf("\nTesting the RK4 and RKF algorythm.\n");
  
  ic = (double *) malloc(sizeof(double) * 4);
  
  ic[0] = 0;
  ic[1] = 100;
  ic[2] = 0;
  ic[3] = 10;

  /* it return the length of the answer vector */
  //n = NUM_rk4 (function, 4, 0.1, 10, ic, &ans, NULL);

  n = NUM_rkf (function, 4, 0.1, 20, ic, &ans, 1e-4, NULL);
  
  for( i = 0; i < n; i++)
  {
    printf("%f  %f  %f  %f  %f\n", ans[4 + 5*i], ans[5*i],
           ans[1+5*i], ans[2 + 5*i], ans[3 + 5*i]);  
  }

  free(ic);
  free(ans);

  return 0;
}

double g1(double x)
{
  return x + 1 - log(x);
}

double f1(double x)
{
  return log(x) - 1;
}

double df1(double x)
{
  return 1/x;
}

double f(double x)
{
  return pow(x-1, 3);
}

int function(int neq, double time, double *y, double *dy, 
             void *data)
{
  
  dy[0] = y[1];
  dy[1] = -9.8;
  dy[2] = y[3];
  dy[3] = 0;
/*
  if (time < 10)
    dy[0] = 9.8 - (13.0/70.0)*y[0];
  else
    dy[0] = 9.8 - (50.0/70.0)*y[0];
*/
  return 0;
}

