#ifndef derivative_h
#define derivative_h

#include "equilibrium.h"

int derivative(equilibrium_t *e);

int fill_temperature_derivative_matrix(double **matrix, equilibrium_t *e);

int fill_pressure_derivative_matrix(double **matrix, equilibrium_t *e);


#endif
