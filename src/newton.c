#include <math.h>
#include <stdio.h>

double NUM_newton(double (*f)(double x), double (*df)(double x), double x0,
                  int nmax, double epsilon)
{

  int i = 0;
  double x1;
  double delta;
  
  do
  {
    if (i >= nmax)
    {
      printf("No convergence within %d iterations.\n", i);
      return 0;
    }
    
    x1 = x0 - f(x0)/df(x0);
    delta = fabs(x1 - x0)/fabs(x1);
    
    x0 = x1;
    i++;
 
  }  while (delta > epsilon);
  printf("Number of iteration: %d\n", i);
  return x1;
}
