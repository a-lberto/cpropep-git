
#ifndef performance_h
#define performance_h

#include "compat.h"
#include "return.h"
#include "equilibrium.h"
#include "derivative.h"



typedef struct _frozen_state_t
{
  float temperature;
  float pressure;
  float velocity;

  float aera_dotm; /* aera / mass flow rate */
  
} frozen_state_t;

typedef struct _equilibrium_state_t
{
  float   temperature;
  float   pressure;
  float   velocity;

  float   aera_dotm;
  
  float   molar_mass;
  deriv_t deriv;
  
} equilibrium_state_t;

typedef struct _performance_t
{
  bool frozen_ok;
  
  struct
  {
    float specific_impulse;
    float molar_mass;

    float cp;     /* specific heat */
    float cp_cv;  /* ratio of specific heat */
    //float isex;

    frozen_state_t chamber;
    frozen_state_t throat;
    frozen_state_t exit;
    
  } frozen;

  bool equilibrium_ok;
  
  struct
  {
    float specific_impulse;

    equilibrium_state_t chamber;
    equilibrium_state_t throat;
    equilibrium_state_t exit;
    
  } equilibrium;

} performance_t;




int frozen_performance(equilibrium_t *e, performance_t *p,
                       double exit_pressure);

int equilibrium_performance(equilibrium_t *e, equilibrium_t *ne,
                            performance_t *p,
                            double exit_pressure);

double mixture_specific_heat_0(equilibrium_t *e, double temp);

#endif

