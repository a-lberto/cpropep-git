/* cpropep.c  -  Calculation of Complex Chemical Equilibrium           */
/* $Id: cpropep.c,v 1.16 2000/05/10 01:35:59 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <ida.pinese@bushnet.qld.edu.au>                      */
/*                                                                     */
/* Licensed under the GPLv2                                            */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>

#include "print.h"
#include "equilibrium.h"
#include "load.h"
#include "libnum.h"

#include "performance.h"

#include "getopt.h"

#include "compat.h"
#include "return.h"

#define version "1.0"
#define date    "17/04/2000"

#define CHAMBER_MSG "Time spent for computing chamber equilibrium"
#define FROZEN_MSG  "Time spent for computing frozen performance"
#define EQUILIBRIUM_MSG "Time spent for computing equilibrium performance"

#define TIME(function, msg) timer = clock(); function;\
                            fprintf(outputfile,\
                                    "%s: %f s\n", msg,\
                                    (double)(clock() - timer)/CLOCKS_PER_SEC)

/* global variable containing the information about chemical species */
extern propellant_t  *propellant_list;
extern thermo_t	     *thermo_list;

extern double conc_tol;
extern double conv_tol;
extern int    iteration_max;

extern FILE * errorfile;
extern FILE * outputfile;

typedef enum _p
{
  SIMPLE_EQUILIBRIUM,
  FIND_FLAME_TEMPERATURE,
  FROZEN_PERFORMANCE,
  EQUILIBRIUM_PERFORMANCE,
  ALL_PERFORMANCE
} p_type;

void welcome_message(void)
{
  printf("----------------------------------------------------------\n");
  printf("Cpropep is an implementation in standard C of the chemical\n"); 
  printf("equilibrium algorythm presented by GORDON and McBRIDE in the\n");
  printf("NASA report SP-273.\n");
  printf("This is the version %s %s\n", version, date);
  printf("This software is release under the GPL and is free of charge\n");
  printf("Copyright (C) 2000 Antoine Lefebvre <antoine.lefebvre@polymtl.ca>\n");
  printf("----------------------------------------------------------\n");
}

void info(char **argv)
{
  printf("Try `%s -h' for more information\n", argv[0]);
}

void usage(void)
{
  printf("Program arguments:\n");
  printf("-f file \t Perform an analysis of the propellant data in file\n");
  printf("-v num  \t Verbosity setting, 0 - 10\n");
  printf("-o file \t Name of the results file, stdout if ommit\n");
  printf("-e file \t Name of the file to store error messages,\
stdout if ommit\n");
  printf("-p      \t Print the propellant list\n");
  printf("-t      \t Print the combustion product list\n");
  printf("-h      \t Print help\n");
  printf("-i      \t Print informations\n");
}


int load_input(FILE *fd, equilibrium_t *e, p_type *p, double *pe)
{ 
  double tmp1, tmp2, m;
  
  int sp;
  int section = 0;

  char buffer[128], tmp[10], num[10], qt[10];
  
  char *bufptr;

  char unit; /* unit of measure of the quantity */

  while ( fgets(buffer, 128, fd) != NULL )
  {
    switch (section)
    {
      case 0:
          
          if (buffer[0] == ' ' || buffer[0] == '\n' || buffer[0] == '\0')
          {
            section = 0;
            break;
          }
          else if ( strncmp(buffer, "Propellant", 10) == 0 )
            section = 1;
          else
          {
            if (strncmp(buffer, "TP", 2) == 0)
              *p = SIMPLE_EQUILIBRIUM;
            else if (strncmp(buffer, "HP", 2) == 0)
              *p = FIND_FLAME_TEMPERATURE;
            else if (strncmp(buffer, "FR", 2) == 0)
              *p = FROZEN_PERFORMANCE;
            else if (strncmp(buffer, "EQ", 2) == 0)
              *p = EQUILIBRIUM_PERFORMANCE;
            else if (strncmp(buffer, "PE", 2) == 0)
              *p = ALL_PERFORMANCE;
            else
            {
              printf ("Unknown option.\n");
              break;
            }
            
            sscanf(buffer, "%s %lf %lf %lf", tmp, &tmp1, &tmp2, pe);
            set_state(e, tmp1, tmp2);
          }
          
          break;
          
      case 1:   /* propellant section */
          if (buffer[0] == '+')
          {
            sscanf(buffer, "%s %s", num, qt);
            bufptr = num + 1;
            sp = atoi(bufptr);
            unit = *(qt + strlen(qt) - 1);
            
            *(qt + strlen(qt) - 1) = '\n';
            
            m = atof(qt);
            
            if ( unit == 'g') 
              add_in_propellant(e, sp, GRAM_TO_MOL(m, sp) );
            else if (unit == 'm')
              add_in_propellant(e, sp, m);
            else
            {
              printf("Unit must be g (gram) or m (mol)\n");
              break;
            }
            break;
          }
          else if (buffer[0] == '#')
          {
            break;
          }
          else if (buffer[0] == ' ' || buffer[0] == '\n' || buffer[0] == '\0')
            section = 0;
          break;
          
      default:
          section = 0;
          break;
    }
  }
  return 0;
}


int main(int argc, char *argv[])
{
  int c, v = 0;
  char filename[FILENAME_MAX];
  FILE *fd = NULL;
  equilibrium_t *equil;

  /* hold performance informations */
  performance_t performance;
  
  
  int thermo_loaded     = 0;
  int propellant_loaded = 0;

  double exit_pressure;

  int temp;

  clock_t timer;
  
  p_type p;

  performance.frozen_ok = false;
  performance.equilibrium_ok = false;
  
  errorfile = stderr;
  outputfile = stdout;
  
//  global_verbose = 1;

  if (argc == 1)
  {
    usage ();
    exit (ERROR);
  }

  
  while (1)
  {
    c = getopt(argc, argv, "ipht?f:v:o:e:");

    if (c == EOF)
      break;
    
    switch (c)
    {
      /* the output file */
      case 'o':
          if (strlen(optarg) > FILENAME_MAX)
          {
            printf("Filename too long!\n");
            break;
          }
          strncpy (filename, optarg, FILENAME_MAX);

          if ( (outputfile = fopen (filename, "w")) == NULL )
            return 1;
          
          break;

          /* the output error file */
      case 'e':
          
          if (strlen(optarg) > FILENAME_MAX)
          {
            printf("Filename too long!\n");
            break;
          }
          strncpy (filename, optarg, FILENAME_MAX);

          if ( (errorfile = fopen (filename, "w")) == NULL )
            return (ERROR);
          
          break;

          /* print the propellant list */
      case 'p':
          if (!propellant_loaded)
          {
            load_propellant ("propellant.dat");
            propellant_loaded = 1;
          }
          print_propellant_list();
          return (SUCCESS);

          /* print the usage */
      case 'h':
          usage();
          return (SUCCESS);

          /* the input file */
      case 'f':
          if (strlen(optarg) > FILENAME_MAX)
          {
            printf("Filename too long!\n");
            break;
          }
          strncpy (filename, optarg, FILENAME_MAX);

          if ( (fd = fopen (filename, "r")) == NULL )
            return (ERROR);
          
          break;

          /* print the thermo list */
      case 't':
          if (!thermo_loaded)
          {
            load_thermo ("thermo.dat");
            thermo_loaded = 1;
          }
          print_thermo_list();
          return (SUCCESS);

          /* set the verbosity level */
      case 'v':
          v = atoi(optarg);
          if (v < 0 || v > 10)
          {
            printf("Verbose is an integer betwenn 0 and 10.\n");
            v = 0;
          }
          break;
          
      case '?':
          info(argv);
          return (SUCCESS);

          /* print information */
      case 'i':
          welcome_message();
          return (SUCCESS);
    }
  }  
  
  if (!thermo_loaded)
  {
    load_thermo ("thermo.dat");
    thermo_loaded = 1;
  }
  if (!propellant_loaded)
  {
    load_propellant ("propellant.dat");
    propellant_loaded = 1;
  }
  
  //welcome_message();


  /* test code to search by formula name
  if ( (temp = propellant_search_by_formula("C")) != -1)
    printf("%s\n", (propellant_list + temp)->name);
  else
    printf("Not found in the database\n");
  */
  
  if (fd != NULL)
  {
    equil = (equilibrium_t *) malloc ( sizeof (equilibrium_t) );
    initialize_equilibrium(equil);
    load_input(fd, equil, &p, &exit_pressure);
    fclose(fd);
    set_verbose(equil, v);

    switch (p)
    {
      case SIMPLE_EQUILIBRIUM:
          print_propellant_composition(equil);
          TIME(if (equilibrium(equil, TP)) break,
               CHAMBER_MSG);
          print_product_composition(equil);
          print_product_properties(equil);
          break;
      case FIND_FLAME_TEMPERATURE:
          print_propellant_composition(equil);
          TIME(if (equilibrium(equil, HP)) break,
               CHAMBER_MSG);
          print_product_composition(equil);
          print_product_properties(equil);
          break;
      case FROZEN_PERFORMANCE:
          print_propellant_composition(equil);
          TIME(if (equilibrium(equil, HP)) break,
               CHAMBER_MSG);
          TIME(frozen_performance(equil, &performance, exit_pressure),
               FROZEN_MSG);
          print_performance_information(&performance);
          break;
      case EQUILIBRIUM_PERFORMANCE:
          print_propellant_composition(equil);
          TIME(if (equilibrium(equil, HP)) break,
               CHAMBER_MSG);
          TIME(equilibrium_performance(equil, &performance, exit_pressure),
               EQUILIBRIUM_MSG);
          print_performance_information(&performance);
          break;
      case ALL_PERFORMANCE:
          print_propellant_composition(equil);
          TIME(if (equilibrium(equil, HP)) break,
               CHAMBER_MSG);
          TIME(frozen_performance(equil, &performance, exit_pressure),
               FROZEN_MSG);
          TIME(equilibrium_performance(equil, &performance, exit_pressure),
               EQUILIBRIUM_MSG);
          print_performance_information(&performance);
          break;
    }
    dealloc_equilibrium (equil);
    free (equil);
  }
  
  free (propellant_list);
  free (thermo_list);

  if (errorfile != stderr)
    fclose (errorfile);

  if (outputfile != stdout)
    fclose (outputfile);
  
  return 0;

}
