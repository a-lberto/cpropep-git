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
#include "libnum.h"

/*from LAPACK */
extern int dgesv_ (int*, int*, double*, int*, int*, double*, int*, int*);
    
int matsol(double **matrix, double *solution, int neq)
{
  int i, j;
  int info;
  
  int    *ipvt; 
  double *mat;
  
  mat  = (double *)malloc(sizeof(double)*neq*neq);
  ipvt = (int *)   malloc(sizeof(int)*neq);

  /* transform the matrix to fit with dgesv_ */
  for (i = 0; i < neq; i++)     
    for(j = 0; j < neq; j++)
      mat[j+neq*i] = matrix[j][i];           

  //copy the right hand side into the vector solution
  for (i = 0; i < neq; i++) 
    solution[i] = matrix[i][neq]; 

  i = 1;
  dgesv_(&neq, &i, mat, &neq, ipvt, solution, &neq, &info);

  //printf("Dgesv return value: %d\n", info);
  /* singular ? */
  if (info > 0)
    return -1;
 
  
  free(ipvt); 
  free(mat);
  
  return 0; 
} 



int gauss(double **matrix, double *solution, int neq)
{   
  int i = 0;
  int j = 0;
  int k = 0;
  double m; /* multiplier */
  double s = 0.0; /* need to store a sommation */
  double temp;
  
  for (k = 0; k < neq-1; k++)
  {
    /* find the smallest j>=k such that a[j][k] != 0 */
    for (j = k; j < neq; j++)
    {
      if (!(matrix[j][k] == 0))
      {
	if (!(j == k))  /* exchange the contents of row j and k */
	{
	  for (i = 0; i <= neq; i++)
	  {
	    temp = matrix[j][i];
	    matrix[j][i] = matrix[k][i];
	    matrix[k][i] = temp;
	  }
	}
	/* the matrix is ready for the next step */
	break;
      }
      printf("No unique solution exists.\n");
      return -1;
    }
    
    for (j = k+1; j < neq; j++)
    {
      /* multiplier */
      m = matrix[j][k]/matrix[k][k];
      
      for (i = k; i <= neq; i++)
      {
        matrix[j][i] = matrix[j][i] - m*matrix[k][i];
      }
    }
  }
  
  if (matrix[neq-1][neq-1] == 0) 
  {
    printf("No unique solution exists.\n");
    return -1;
  }
  else
  {
    solution[neq-1] = matrix[neq-1][neq]/matrix[neq-1][neq-1];
    
    for (j = neq-2; j >= 0; j--)
    {
      s = 0.0;
      for (i = j + 1; i < neq; i++)
        s += matrix[j][i]*solution[i];
      
      solution[j] = (matrix[j][neq] - s)/matrix[j][j];
    }
  }
  return 0;
}

/* LU-Factorisation, Doolittle's method */
int lu(double **matrix, double *solution, int neq)
{
  int i, j, s;
  double tmp = 0.0;
  
  /* L is a lower triangular matrix with diagonal set to one */
  double **L;
  /* U is an upper triangular matrix */
  double **U;
  
  /* temporary vector Ly = b, Ux = y */
  double *y;
  
  
  y = (double *)calloc(neq, sizeof(double));
  
  L = (double **)calloc(neq, sizeof(double*));
  U = (double **)calloc(neq, sizeof(double*));
  
  
  for (i = 0; i < neq; i++)
  {
    L[i] = (double *)calloc(neq, sizeof(double));
    U[i] = (double *)calloc(neq, sizeof(double));
    solution[i] = 0; /* reset the solution vector */
  }
  
  /* set the diagonal to 1 */
  for (i = 0; i < neq; i++)
    L[i][i] = 1;

  
  /* LU Factorisation */
  for (i = 0; i < neq; i++)
  {
    U[0][i] = matrix[0][i];
	 
    if (i > 0)
    {
      for (j = 1; j <= i; j++)
      {
	tmp = 0.0;
	for (s = 0; s < j; s++)
	  tmp += L[j][s]*U[s][i];
	
	U[j][i] = matrix[j][i] - tmp;
      }
    }
	 
    for (j = i + 1; j < neq; j++)
    {
      if (U[i][i] == 0.0)
      {
	printf("No unique solution exist.\n");
	return -1;
      }
      if (i == 0)
	L[j][i] = matrix[j][i]/U[i][i];
      
      else
      {
	tmp = 0.0;
	for (s = 0; s < i; s++)
	  tmp += L[j][s]*U[s][i];
	
	L[j][i] = (matrix[j][i] - tmp)/U[i][i];
      }
    }
  }
  /* End LU-Factorisation */
  
  /* substitution  for y    Ly = b*/
  for (i = 0; i < neq; i++)
  {
    tmp = 0.0;
    
    for (j = 0; j < i; j++)
      tmp += L[i][j]*y[j];
    
    y[i] = matrix[i][neq] - tmp;
  }
  
  /* substitution for x   Ux = y*/
  for (i = neq-1; i >=0; i--)
  {
    if (U[i][i] == 0.0)
    {
      printf("No unique solution exist.\n");
      return -1;
    }
    
    tmp = 0.0;
    for (j = i; j < neq; j++)
      tmp += U[i][j]*solution[j];
    
    solution[i] = (y[i] - tmp)/U[i][i];
  }
     
  /*print_square_matrix(L, neq);*/
  /*print_square_matrix(L, neq);*/
     
  free (y);
  free (L);
  free (U);
  return 0;      
}

int print_square_matrix(double **matrix, int neq)
{
  int i = 0;
  int j = 0;
  
  for (i = 0; i < neq; i++)
  {
    for (j = 0; j < neq; j++)
      printf("% .5f ", matrix[i][j]);
    printf("\n");
  }
  printf("\n");
  return 0;
}

int print_matrix(double **matrix, int neq)
{
  int i = 0;
  int j = 0;
  
  for (i = 0; i < neq; i++)
  {
    for (j = 0; j <= neq; j++)
      printf("% .5e ", matrix[i][j]);
    printf("\n");
  }
  printf("\n");
  return 0;
}

int print_vec(double *vec, int neq)
{
  int i;
  for (i = 0; i < neq; i++)
    printf("% .5e ", vec[i]);
  printf("\n");
  return 0;
}

int rk4( int (*f)(int neq, double time, double *y, double *dy, 
		  void *data), 
	 int neq, double step, double duration, double *ic, 
	 double **y, void *data)
{
  int i;
  int n;

  double t = 0.0;

  double *tmp;
  double *dy;
  double *K1, *K2, *K3, *K4;   

  tmp = (double *)malloc(sizeof(double) * neq);
  dy = (double *)malloc(sizeof(double) * neq);

  K1 = (double *)malloc(sizeof(double) * neq);
  K2 = (double *)malloc(sizeof(double) * neq);
  K3 = (double *)malloc(sizeof(double) * neq);
  K4 = (double *)malloc(sizeof(double) * neq);

  for (i = 0; i < neq; i++)
  {
    y[0][i] = ic[i]; // conditions initiales
    tmp[i] = y[0][i];
  }
 
  for (n = 0; n < (int)round(duration/step); n++)
  {

    for (i = 0; i < neq; i++)
    {
      f(neq, t, tmp, dy, data);
      K1[i] = step*dy[i];
      
      tmp[i] = y[n][i] + K1[i]/2;  // for the next step           
    }
    
    for (i = 0; i < neq; i++)
    {
      f(neq, t, tmp, dy, data);
      K2[i] = step*dy[i];
      
      tmp[i] = y[n][i] + K2[i]/2;
    }
    
    for (i = 0; i < neq; i++)
    {
      f(neq, t, tmp, dy, data);
      K3[i] = step*dy[i];
      
      tmp[i] = y[n][i] + K3[i];
    }
    
    for (i = 0; i < neq; i++)
    {
      f(neq, t, tmp, dy, data);
      K4[i] = step*dy[i];
    }
    
    for (i = 0; i < neq; i++)
      y[n+1][i] = y[n][i] + (1.0/6.0)*(K1[i] + 2.0*K2[i] + 2.0*K3[i] + K4[i]);

    t = t + step;
  }

  free(tmp);
  free(dy);
  free(K1);
  free(K2);
  free(K3);
  free(K4);
  return 0;
}



int round(double a)
{
  int t = a;
  
  if (a - (double)t < 0.5)
    return t;
  else if (a - (double)t > 0.5)
    return t + 1;
  else
    return (t + (t % 2));
}




