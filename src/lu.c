#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "num.h"

/* LU-Factorisation, Doolittle's method */
int NUM_lu(double *matrix, double *solution, int neq)
{
  int i, j, s;
  double tmp = 0.0;
  
  /* L is a lower triangular matrix with diagonal set to one */
  double *L;
  /* U is an upper triangular matrix */
  double *U;  
  /* temporary vector Ly = b, Ux = y */
  double *y;
  
  y = (double *) calloc (neq, sizeof(double));
  L = (double *) calloc (neq*neq, sizeof(double));
  U = (double *) calloc (neq*neq, sizeof(double));
  
  for (i = 0; i < neq; i++)
  {
    solution[i]  = 0; /* reset the solution vector */
    L[i + neq*i] = 1;
  }
    
  /* LU Factorisation */
  for (i = 0; i < neq; i++)
  {
    U[0 + neq*i] = matrix[0 + neq*i];
    
    if (i > 0)
    {
      for (j = 1; j <= i; j++)
      {
        tmp = 0.0;
        for (s = 0; s < j; s++)
          tmp += L[j + neq*s] * U[s + neq*i];
        
        U[j + neq*i] = matrix[j + neq*i] - tmp;
      }
    }
    
    for (j = i + 1; j < neq; j++)
    {
      if (U[i + neq*i] == 0.0)
      {
        printf("LIBNUM: No unique solution exist.\n");
        return -1;
      }
      if (i == 0)
      {
        L[j + neq*i] = matrix[j + neq*i]/U[i + neq*i];      
      }
      else
      {
        tmp = 0.0;
        for (s = 0; s < i; s++)
          tmp += L[j + neq*s]*U[s + neq*i];
        
        L[j + neq*i] = (matrix[j + neq*i] - tmp)/U[i + neq*i];    
      }
    }
  }
  /* End LU-Factorisation */
  
  /* substitution  for y    Ly = b*/
  for (i = 0; i < neq; i++)
  {
    tmp = 0.0;
    for (j = 0; j < i; j++)
      tmp += L[i + neq*j]*y[j];
    
    y[i] = matrix[i + neq*neq] - tmp;
  }
  
  /* substitution for x   Ux = y*/
  for (i = neq - 1; i >=0; i--)
  {
    if (U[i + neq*i] == 0.0)
    {
      printf("LIBNUM: No unique solution exist.\n");
      return -1;
    }
    
    tmp = 0.0;
    for (j = i; j < neq; j++)
      tmp += U[i + neq*j]*solution[j];
    
    solution[i] = (y[i] - tmp)/U[i + neq*i];

    /* isnan() is probably not available under Windows */
#ifdef LINUX
    if ( isnan(solution[i]) )
    {
      printf("LIBNUM: No unique solution exist. NaN in solution.\n");
      return -1;
    }
#endif
    
  }
     
  free (y);
  free (L);
  free (U);
  return 0;      
}

/* This function will divide each row of the matrix
 * by the highest element of this row.
 * This kind of scaling could improve precision while
 * solving some difficult matrix.
 */
int NUM_matscale(double *matrix, int neq)
{
  int i; /* line */
  int j; /* column */

  double val;
  double tmp;
  
  for (i = 0; i < neq; i++)
  {
    val = 0;
    /* find the highest value */
    for (j = 0; j < neq; j++)
    {
      tmp = abs(matrix[i + neq*j]);
      val = (tmp > val) ? tmp : val;
    }

    /* divide element of the line by this value
     * including the right side  */
    if (val == 0)
      return -1;
    
    for (j = 0; j < neq+1; j++)
      matrix[i + neq*j] = matrix[i + neq*j]/val;
  }
  return 0;

}
