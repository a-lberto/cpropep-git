
#ifndef performance_h
#define performance_h

#include "equilibrium.h"

int frozen_performance(equilibrium_t *e, double exit_pressure);
int equilibrium_performance(equilibrium_t *e, double exit_pressure);

double mixture_specific_heat_0(equilibrium_t *e, double temp);

#endif
