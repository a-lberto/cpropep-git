#ifndef derivative_h
#define derivative_h

#include "equilibrium.h"
typedef struct deriv
{
  double del_lnV_lnT; /* derivative of ln(V) with respect to ln(T)
                         at constant pressure */
  double del_lnV_lnP; /* derivative of ln(V) with respect to ln(P)
                         at constant temperature */
  double cp;
  double cv;
  double cp_cv;
  double isex;        /* isentropic exponent */
  
} deriv_t;

int derivative(equilibrium_t *e, deriv_t *d);

int fill_temperature_derivative_matrix(double **matrix, equilibrium_t *e);

int fill_pressure_derivative_matrix(double **matrix, equilibrium_t *e);


#endif
