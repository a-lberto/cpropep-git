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

int main(void)
{

  /*
  double *matrix;
  double *solution;
  
  matrix = (double *) malloc (sizeof(double)*12);
  solution = (double *) malloc (sizeof(double)*3);
  
  matrix[0] = 1;
  matrix[1] = 1;
  matrix[2] = 9;
  matrix[3] = -1;
  matrix[4] = 1;
  matrix[5] = 3;
  matrix[6] = 1;
  matrix[7] = 1;
  matrix[8] = 1;
  matrix[9] = 6;
  matrix[10] = 2;
  matrix[11] = 22;

  lu(matrix, solution, 3);

  print_vec(solution, 3);
  */

  printf("Solution: %f\n",  sec(f, -15, 0, 0.0001));
    
  //printf("epsilon: %.16e\n", epsilon());
  
  /*
  int i;
  double **ans;
  double *ic;
  
  ic = (double *)malloc(sizeof(double) * 4);
  
  ic[0] = 0;
  ic[1] = 100;
  ic[2] = 0;
  ic[3] = 10;

  if ( ( ans = (double **)malloc(sizeof(double *) * 101)) == NULL)
       printf("Problem allocating memory\n");
  for (i = 0; i < 101; i++)
    ans[i] = (double *)malloc(sizeof(double) * 4);
 
  rk4 (function, 4, 0.1, 10, ic, ans, NULL);

  for( i = 0; i < 100; i++)
    printf("%f  %f %f %f %f\n", i*0.1, ans[i][0], ans[i][1], 
	   ans[i][2], ans[i][3]);  

  return 0;


  free(ic);
  free(ans);
  */

  return 0;
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

  return 0;
}

