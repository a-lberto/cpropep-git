
#include <math.h>
#include <stdio.h>

double ptfix(double (*f)(double x), double x0, double nmax, double epsilon)
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
    
    x1 = f(x0);

    delta = fabs(x1 - x0)/fabs(x1);
    x0 = x1;
    i++;

  } while (delta > epsilon);
  printf("Number of iteration: %d\n", i);
  return x1;
}
