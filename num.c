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

//#include <time.h>

#include "libnum.h"

int function(int neq, double time, double *y, double *dy, 
	     int ierr);

int main(void)
{

  /* for timing */
  /* clock_t start, finish; */

  
  //const int neq = 3;
  //int i = 0;
  
  //double **mat;
  //double *ans;


  //mat = (double **)calloc(neq, sizeof(double*));
  //for (i = 0; i < neq; i++)
  //{
  //  mat[i] = (double *)calloc(neq+1, sizeof(double));
  //}
  //ans = (double *)calloc(neq, sizeof(double));

  //mat[0][0] = 4;
  //mat[0][1] = 10;
  //mat[0][2] = -2;
  //mat[0][3] = -20;
  /* mat[0][4] = 54; */

  //mat[1][0] = -1;
  //mat[1][1] = -15;
  //mat[1][2] = 3;
  //mat[1][3] = 30;
  /* mat[1][4] = 0; */
  
  //mat[2][0] = 0;
  //mat[2][1] = 25;
  //mat[2][2] = -5;
  //mat[2][3] = -50;
  /* mat[2][4] = 2; */

  /* mat[3][0] = 43; */
  /* mat[3][1] = 0; */
  /* mat[3][2] = 76; */
  /* mat[3][3] = 0; */
  /* mat[3][4] = 77; */

  //print_matrix(mat, neq); 
  
  //lu(mat, ans, neq);
  
  /* gauss(mat, ans, neq); */
     
  //print_matrix(mat, neq);
  //print_vec(ans, neq);
  
  //free(mat);
  //free(ans);

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
 
  

  rk4 (function, 4, 0.1, 10, ic, ans);

  for( i = 0; i < 100; i++)
    printf("%f  %f %f %f %f\n", i*0.1, ans[i][0], ans[i][1], 
	   ans[i][2], ans[i][3]);  

  return 0;


  free(ic);
  free(ans);

}


//double ffunction(double t, double y)
//{
//  return ( (t - y) / 2.0 );
//}

int function(int neq, double time, double *y, double *dy, 
	     int ierr)
{
  
  dy[0] = y[1];
  dy[1] = -9.8;
  dy[2] = y[3];
  dy[3] = 0;

  return 0;
}

