#include <stdlib.h>
#include <math.h>
#include "num.h"

/* neq:      number of equations
 * step:     initial time step
 * duration: total duration of the simulation
 * ic:       initial conditions vector
 * y:        solution vector
 * epsil:    precision on the data
 * data:     pointer to data structure
 */

/* Runge-Kutta-Fehlberg 4-5 order */
 
int rkf(int (*f)(int neq, double time, double *y, double *dy, void *data), 
        int neq, double step, double duration, double *ic, 
        double *y, double epsil, void *data)
{
  int i;
  int n;

  double h;
  double t = 0.0;

  double *tmp;
  double *dy;
  double *K1, *K2, *K3, *K4, *K5, *K6;   

  double *E; /* Error vector of the error */
  
  double err; /* maximum error */
  double beta;
  
  tmp = (double *) malloc(sizeof(double) * neq);
  dy = (double *) malloc(sizeof(double) * neq);

  E = (double *) malloc(sizeof(double) * neq);
  
  K1 = (double *) malloc(sizeof(double) * neq);
  K2 = (double *) malloc(sizeof(double) * neq);
  K3 = (double *) malloc(sizeof(double) * neq);
  K4 = (double *) malloc(sizeof(double) * neq);
  K5 = (double *) malloc(sizeof(double) * neq);
  K6 = (double *) malloc(sizeof(double) * neq);

  h = step;
  n = 0;
  y = (double *) malloc(sizeof(double) * neq * n);
  
  for (i = 0; i < neq; i++)
  {
    y[i + neq*0] = ic[i]; // conditions initiales
    tmp[i] = y[i + neq*0];
  }

  while (t < duration)
  {

    y = (double *) realloc(y, sizeof(double) * neq * (n + 1));
    

    for (i = 0; i < neq; i++)
    {
      f(neq, t, tmp, dy, data);
      K1[i] = h * dy[i];
      
      tmp[i] = y[i + neq*n] + K1[i]/4;  // for the next step           
    }
      
    for (i = 0; i < neq; i++)
    {
      f(neq, t + h/4, tmp, dy, data);
      K2[i] = h * dy[i];

      tmp[i] = y[i + neq*n] + 3*K1[i]/32 + 9*K2[i]/32;
        
    }

    for (i = 0; i < neq; i++)
    {
      f(neq, t + 3*h/8, tmp, dy, data);
      K3[i] = h * dy[i];

      tmp[i] = y[i + neq*n] + 1932*K1[i]/2197 - 7200*K2[i]/2197 +
        7296*K3[i]/2197;
    }

    
    for (i = 0; i < neq; i++)
    {
      f(neq, t + 12*h/13, tmp, dy, data);
      K4[i] = h * dy[i];
      
      tmp[i] = y[i + neq*n] + 439*K1[i]/216 - 8*K2[i] + 3680*K3[i]/513 -
        845*K4[i]/4104;
    }

    for (i = 0; i < neq; i++)
    {
      f(neq, t + h, tmp, dy, data);
      K5[i] = h * dy[i];
      
      tmp[i] = y[i + neq*n] + 8*K1[i]/27 + 2*K2[i] - 3544*K3[i]/2565 +
        1859*K4[i]/4104 - 11*K5[i]/40;
    }
    
    for (i = 0; i < neq; i++)
    {
      f(neq, t + h/2, tmp, dy, data);
      K6[i] = h * dy[i];
    }

    err = 0;
    for (i = 0; i <  neq; i++)
    {
      E[i] = fabs((K1[i]/360 - 128*K3[i]/4275 - 2197*K4[i]/75240 +
                   K5[i]/50 + 2*K6[i]/55)/h);
      err = ((E[i] > err) ? E[i] : err);      
    }

    if (err < epsil)
    {

      for (i = 0; i < neq; i++)
      {
        y[i + neq*(n+1)] = y[i + neq*n] + 25*K1[i]/216 + 1408*K3[i]/2565 +
          2197*K4[i]/4104 - 0.2*K5[i];
      }
      t += h;
      n++;
    }
    else
    {
      beta = pow(epsil/(2*err), 0.25);

      if (beta < 0.1)
      {
        h = 0.1* h;
      }
      else if (beta > 4)
      {
        h = 4*h;
      }
      else
      {
        h = beta * h;
      }
    }

    if (t + h > duration)
    {
      h =  duration - t;
    }   
  }

  return 0;
}

