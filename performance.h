
#ifndef performance_h
#define performance_h

#include "equilibrium.h"

#include "compat.h"

typedef struct _performance_t
{

  bool frozen_ok;
  
  struct
  {
    float specific_impulse;

    float cp;     /* specific heat */
    float cp_cv;  /* ratio of specific heat */
    //float isex;
    float molar_mass;
    
    struct 
    {
      float temperature;
      float pressure;
      float velocity;
    } throat;

    struct
    {
      float temperature;
      float pressure;
      float velocity;
    } exit;

  } frozen;

  bool equilibrium_ok;
  
  struct
  {
    float specific_impulse;

    //float cp;     /* specific heat */
    //float cp_cv;  /* ratio of specific heat */
    
    struct 
    {
      float temperature;
      float pressure;
      float velocity;

      float cp;
      float cp_cv;
      float isex;
      float molar_mass;
      
    } throat;

    struct
    {
      float temperature;
      float pressure;
      float velocity;

      float cp;
      float cp_cv;
      float isex;
      float molar_mass;
      
    } exit;
    
  } equilibrium;

} performance_t;




int frozen_performance(equilibrium_t *e, performance_t *p,
                       double exit_pressure);

int equilibrium_performance(equilibrium_t *e, performance_t *p,
                            double exit_pressure);

double mixture_specific_heat_0(equilibrium_t *e, double temp);

#endif
