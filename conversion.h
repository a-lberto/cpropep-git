#ifndef conversion_h
#define conversion_h


/* pressure units */
enum
{
  ATM,
  PSI,
  BAR,
  KPA
};

/* Transform calories to joules */
#define CAL_TO_JOULE       4.1868

/* Transform pound/(cubic inch) to gram/(cubic centimeter) */
#define LBS_IN3_TO_G_CM3  27.679905

#define ATM_TO_PA         101325
#define ATM_TO_PSI        14.695949
#define ATM_TO_BAR         1.01325

#define BAR_TO_PSI        14.503774

#define BAR_TO_ATM         0.98692327
#define PSI_TO_ATM         0.068045964
#define KPA_TO_ATM         0.0098692327

#endif
