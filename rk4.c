#include <stdlib.h>
#include "num.h"

int rk4(int (*f)(int neq, double time, double *y, double *dy, void *data), 
        int neq, double step, double duration, double *ic, 
        double *y, void *data)
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
    y[i + neq*0] = ic[i]; // conditions initiales
    tmp[i] = y[i + neq*0];
  }
 
  for (n = 0; n < Round(duration/step); n++)
  {

    for (i = 0; i < neq; i++)
    {
      f(neq, t, tmp, dy, data);
      K1[i] = step * dy[i];
      
      tmp[i] = y[i + neq*n] + K1[i]/2;  // for the next step           
    }

    /* FIXME: verify the coefficient t + step/2 */
    for (i = 0; i < neq; i++)
    {
      f(neq, t + step/2, tmp, dy, data);
      K2[i] = step * dy[i];
      
      tmp[i] = y[i + neq*n] + K2[i]/2;
    }
    
    for (i = 0; i < neq; i++)
    {
      f(neq, t + step/2, tmp, dy, data);
      K3[i] = step * dy[i];
      
      tmp[i] = y[i + neq*n] + K3[i];
    }
    
    for (i = 0; i < neq; i++)
    {
      f(neq, t + step, tmp, dy, data);
      K4[i] = step * dy[i];
    }
    
    for (i = 0; i < neq; i++)
      y[i + neq*(n+1)] = y[i + neq*n] +
        (1.0/6.0)*(K1[i] + 2.0*K2[i] + 2.0*K3[i] + K4[i]);

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
