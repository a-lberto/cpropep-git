/* test.c - Testing the functionnality of the various method
 * $Id: test.c,v 1.3 2000/10/20 20:17:20 antoine Exp $
 * Copyright (C) 2000
 *    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>
 *
 * Licensed under the GPL
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include <time.h>

#include "num.h"

FILE * errorfile;
FILE * outputfile;

int test_rk4(void);
int test_lu(void);
int test_sysnewton(void);
int test_sec(void);
int test_newton(void);
int test_ptfix(void);

/* g1(x) = x + 1 - ln(x) */
double g1(double x) {
  return x + 1 - log(x); 
}
/* f1(x) = ln(x) - 1 */
double f1(double x) {
  return log(x) - 1; 
}
/* f1'(x) = 1/x */
double df1(double x) {
  return 1/x; 
}
double dg1(double x) {
  return 1 - 1/x;
}
/* f(x) = (x-1)^3 */
double f(double x) {
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

/* functions to test the sysnewton algorythm
 *
 * The system to be solve is the following
 * r1(x1, x2) = e^x1 - x2 = 0
 * r2(x1, x2) = x1^2 + x2^2 - 16 = 0
 *
 */
double r1(double *x)
{
  return exp(x[0]) - x[1];
}
double r2(double *x)
{
  return pow(x[0], 2) + pow(x[1], 2) - 16;
}
/* We need partial derivative of these function with each variable
 * dr1_dx1(x1, x2) = e^x1
 * dr1_dx2(x1, x2) = -1
 * dr2_dx1(x1, x2) = 2*x1
 * dr2_dx2(x1, x2) = 2*x2
 */
double dr1_dx1(double *x)
{
  return exp(x[0]);
}
double dr1_dx2(double *x)
{
  return -1;
}
double dr2_dx1(double *x)
{
  return 2*x[0];
}
double dr2_dx2(double *x)
{
  return 2*x[1];
}


int main(void)
{

  test_lu();
//  test_rk4();
  test_sysnewton();

  test_sec();
  test_newton();
  test_ptfix();
 
  return 0;
}

int test_sec(void)
{
  double ans;
  printf("Testing Secante method\n");
  if (NUM_sec(f1, 1, 2, 100, 0.0001, &ans))
    printf("No solution: error in the method.\n\n");
  else
    printf("Solution: %f\n\n",  ans);
  return 0;
}

int test_newton(void)
{
  double ans;
  printf("Testing Newton method\n");
  if (NUM_newton(f1, df1, 1, 100, 0.0001, &ans))
    printf("No solution: error in the method.\n\n");
  else
    printf("Solution: %f\n\n", ans);
  return 0;
}

int test_ptfix(void)
{
  double ans;
  printf("Testing fixed point method\n");
  if (NUM_ptfix(g1, 1, 100, 0.0001, &ans))
     printf("No solution: error in the method.\n");
  else
    printf("Solution: %f\n\n", ans);
  return 0;
}

int test_sysnewton(void)
{
  func_t *jac;
  func_t *r;
  int nvar = 2;

  double x[2]; /* solution vector, contain initially the initial conditions */
  
  printf("Testing newton method for solving"
         " non-linear system of equations.\n");

  jac = (func_t *) malloc (sizeof(func_t) * nvar * nvar);
  r = (func_t *) malloc (sizeof(func_t) * nvar);

  /* Initialize the functions pointers array */
  r[0] = r1;
  r[1] = r2;

  jac[0] = dr1_dx1;
  jac[1] = dr2_dx1;
  jac[2] = dr1_dx2;
  jac[3] = dr2_dx2;

  /* set the initial estimate */
  x[0] = 2.8;
  x[1] = 2.8;

  /* call the sysnewton function */

  NUM_sysnewton(jac, r, x, 2, 100, 1e-8);

  /* print the solution */

  printf("Solution: x1 = %f, x2 = %f\n", x[0], x[1]);

  printf("\n");
  return 0;
  
}

int test_lu(void)
{
  int i;
  double *matrix;
  double *solution;
  int size = 4;
  
  printf("Testing the LU factorisation algotythm.\n");
  matrix = (double *) malloc (sizeof(double)*size*(size+1));
  solution = (double *) malloc (sizeof(double)*size);
  
  matrix[0] = 1;
  matrix[1] = 1;
  matrix[2] = 1;
  matrix[3] = 1;

  matrix[4] = 1;
  matrix[5] = 2;
  matrix[6] = 1;
  matrix[7] = 1;

  matrix[8] = 1;
  matrix[9] = 1;
  matrix[10] = 0;
  matrix[11] = -1;
  
  matrix[12] = 1;
  matrix[13] = 1;
  matrix[14] = 1;
  matrix[15] = -1;
  
  matrix[16] = 4;
  matrix[17] = 5;
  matrix[18] = 3;
  matrix[19] = 0;

  //NUM_matscale(matrix, size);
  NUM_print_matrix(matrix, size);

  if (NUM_lu(matrix, solution, size))
    printf("No solution: Error in the numerical method,\n");
  else
    NUM_print_vec(solution, size);

  for (i = 0; i < size; i++)
  {
    if (solution[i] != 1.0)
    {
      printf("Error found in the solution.\n");
    }
  }
  printf("\n");
  return 0;
}


int test_rk4(void)
{
  int i, n;
  double *ans;
  double *ic;
  
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
