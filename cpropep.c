/* cpropep.c  -  Calculation of Complex Chemical Equilibrium           */
/* $Id: cpropep.c,v 1.22 2000/07/11 03:48:57 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <pinese@cyberwizards.com.au>                         */
/*                                                                     */
/* Licensed under the GPLv2                                            */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>

#include "getopt.h"
#include "load.h"

#include "equilibrium.h"
#include "performance.h"
#include "derivative.h"
#include "thermo.h"

#include "print.h"

#include "compat.h"
#include "return.h"

#define version "1.0"
#define date    "10/07/2000"

#define CHAMBER_MSG     "Time spent for computing chamber equilibrium"
#define FROZEN_MSG      "Time spent for computing frozen performance"
#define EQUILIBRIUM_MSG "Time spent for computing equilibrium performance"

#define TIME(function, msg) timer = clock(); function;\
                            fprintf(outputfile,\
                                    "%s: %f s\n\n", msg,\
                                    (double)(clock() - timer)/CLOCKS_PER_SEC)

#define THERMO_FILE     "thermo.dat"
#define PROPELLANT_FILE "propellant.dat"

#undef TIME
#define TIME(function, msg) function;

#define MAX_CASE 10

typedef enum _p
{
  SIMPLE_EQUILIBRIUM,
  FIND_FLAME_TEMPERATURE,
  FROZEN_PERFORMANCE,
  EQUILIBRIUM_PERFORMANCE
} p_type;

char case_name[][80] = {
  "Fixed pressure-temperature equilibrium",
  "Fixed enthalpy-pressure equilibrium - adiabatic flame temperature",
  "Frozen equilibrium performance evaluation",
  "Shifting equilibrium performance evaluation"
};

typedef struct _case_t
{
  p_type p;

  bool temperature_set;
  bool pressure_set;
  bool exit_condition_set;

  double           temperature;
  double           pressure;
  exit_condition_t exit_cond_type;
  double           exit_condition;
  
} case_t;


void welcome_message(void)
{
  printf("----------------------------------------------------------\n");
  printf("Cpropep is an implementation in standard C of the chemical\n"); 
  printf("equilibrium algorythm presented by GORDON and McBRIDE in the\n");
  printf("NASA report RP-1311.\n");
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
  printf("-q num  \t Print information about propellant number num\n");
  printf("-t      \t Print the combustion product list\n");
  printf("-u num  \t Print information about product number num\n");
  printf("-h      \t Print help\n");
  printf("-i      \t Print informations\n");
}


int load_input(FILE *fd, equilibrium_t *e, case_t *t, double *pe)
{ 
  double m;
  
  int sp;
  int section = 0;

  int n_case = 0;
  
  char buffer[128], num[10], qt[10], unit[10];
  char variable[64];
  
  char *bufptr;

  while ( fgets(buffer, 128, fd) != NULL )
  {
    switch (section)
    {
      case 0:

          if (n_case >= MAX_CASE)
          {
            fprintf(outputfile,
          "Warning: Too many different case, maximum is %d: deleting case.\n",
                   MAX_CASE+1);
            section = 100;
            break;
          }
          
          if (buffer[0] == ' ' || buffer[0] == '\n' || buffer[0] == '\0' ||
              buffer[0] == '#')
          {
            section = 0;
            break;
          }
          else if ( strncmp(buffer, "Propellant", 10) == 0 )
          {
            section = 1;
          }
          else
          { 
            if (strncmp(buffer, "TP", 2) == 0)
              t[n_case].p = SIMPLE_EQUILIBRIUM;
            else if (strncmp(buffer, "HP", 2) == 0)
              t[n_case].p = FIND_FLAME_TEMPERATURE;
            else if (strncmp(buffer, "FR", 2) == 0)
              t[n_case].p = FROZEN_PERFORMANCE;
            else if (strncmp(buffer, "EQ", 2) == 0)
              t[n_case].p = EQUILIBRIUM_PERFORMANCE;
            else
            {
              printf ("Unknown option.\n");
              break;
            }
            section = 2;
          }
          
          break;
          
      case 1:   /* propellant section */
          if (buffer[0] == '+')
          {
            sscanf(buffer, "%s %s %s", num, qt, unit);
            bufptr = num + 1;
            sp = atoi(bufptr);
            
            m = atof(qt);
            
            if (strcmp(unit, "g") == 0)
            {
              add_in_propellant(e, sp, GRAM_TO_MOL(m, sp) );
            }
            else if (strcmp(unit, "m") == 0)
            {
              add_in_propellant(e, sp, m);
            }
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
          {
            section = 0;
          }
          break;

      case 2:
          if (buffer[0] == '+')
          {
            sscanf(buffer, "%s %s %s", variable, qt, unit);

            bufptr = variable + 1;

            if (strcmp(bufptr, "chamber_temperature") == 0)
            {
              t[n_case].temperature     = atof(qt);
              t[n_case].temperature_set = true;
            }
            else if (strcmp(bufptr, "chamber_pressure") == 0)
            {
              t[n_case].pressure     = atof(qt);
              t[n_case].pressure_set = true;
            }
            else if (strcmp(bufptr, "exit_pressure") == 0)
            {
              t[n_case].exit_cond_type     = PRESSURE;
              t[n_case].exit_condition     = atof(qt);
              t[n_case].exit_condition_set = true;
              
            }
            else if (strcmp(bufptr, "supersonic_area_ratio") == 0)
            {
              t[n_case].exit_cond_type = SUPERSONIC_AREA_RATIO;
              t[n_case].exit_condition = atof(qt);
              t[n_case].exit_condition_set = true;
            }
            else if (strcmp(bufptr, "subsonic_area_ratio") == 0)
            {
              t[n_case].exit_cond_type = SUBSONIC_AREA_RATIO;
              t[n_case].exit_condition = atof(qt);
              t[n_case].exit_condition_set = true;
            }
            else
            {
              printf("Unknown keyword.\n");
              break;
            }
 
            break;

          }
          else if (buffer[0] == '#')
          {
            break;
          }
          else if (buffer[0] == ' ' || buffer[0] == '\n' || buffer[0] == '\0')
          {
            section = 0;
            n_case++;
          }
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
  int i, c, v = 0;
  char filename[FILENAME_MAX];
  FILE *fd = NULL;

  equilibrium_t *equil, *frozen, *shifting; 
  
  int thermo_loaded     = 0;
  int propellant_loaded = 0;

  double exit_pressure;

  clock_t timer;

//  int param;

  case_t case_list[MAX_CASE];
  for (i = 0; i < MAX_CASE; i++)
  {
    case_list[i].p = -1;
    case_list[i].temperature_set = false;
    case_list[i].pressure_set = false;
    case_list[i].exit_condition_set = false;
  }
  
  errorfile = stderr;
  outputfile = stdout;
  
/*  global_verbose = 1; */

  if (argc == 1)
  {
    usage ();
    exit (ERROR);
  }

  
  while (1)
  {
    c = getopt(argc, argv, "ipht?f:v:o:e:q:u:");

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
            load_propellant (PROPELLANT_FILE);
            propellant_loaded = 1;
          }
          print_propellant_list();
          free(propellant_list);
          return (SUCCESS);

          /* print propellant info */
      case 'q':
          if (!propellant_loaded)
          {
            load_propellant (PROPELLANT_FILE);
            propellant_loaded = 1;
          }
          print_propellant_info( atoi(optarg) );
          free(propellant_list);
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
            load_thermo (THERMO_FILE);
            thermo_loaded = 1;
          }
          print_thermo_list();
          free(thermo_list);
          return (SUCCESS);

      case 'u':
          if (!thermo_loaded)
          {
            load_thermo (THERMO_FILE);
            thermo_loaded = 1;
          }
          print_thermo_info( atoi(optarg) );
          free(thermo_list);
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
    load_thermo (THERMO_FILE);
    thermo_loaded = 1;
  }
  if (!propellant_loaded)
  {
    load_propellant (PROPELLANT_FILE);
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
    equil = (equilibrium_t *) malloc (sizeof (equilibrium_t));
    initialize_equilibrium(equil);


    frozen   = (equilibrium_t *) malloc (sizeof(equilibrium_t)*3);
    shifting = (equilibrium_t *) malloc (sizeof(equilibrium_t)*3);

    for (i = 0; i < 3; i++)
    {
      initialize_equilibrium(frozen + i);
      initialize_equilibrium(shifting + i);
    }
    
    
    load_input(fd, equil, case_list, &exit_pressure);

    compute_density(&(equil->propellant));
    
    fclose(fd);
    global_verbose = v;

    list_element(equil);
    list_product(equil);
    
    i = 0;
    while ((case_list[i].p != -1) && (i <= MAX_CASE))
    {
      fprintf(outputfile, "Computing case %d\n%s\n\n", i+1,
              case_name[case_list[i].p]);
      
      switch (case_list[i].p)
      {
        case SIMPLE_EQUILIBRIUM:

            if (!(case_list[i].temperature_set))
            {
              printf("Chamber temperature not set. Aborted.\n");
              break;
            }
            else if (!(case_list[i].pressure_set))
            {
              printf("Chamber pressure not set. Aborted.\n");
              break;
            }

            equil->properties.T = case_list[i].temperature;
            equil->properties.P = case_list[i].pressure;
            
            print_propellant_composition(equil);
            TIME(if (equilibrium(equil, TP)) break,
                 CHAMBER_MSG);

            print_product_properties(equil, 1);
            print_product_composition(equil, 1);
            break;

        case FIND_FLAME_TEMPERATURE:

            if (!(case_list[i].pressure_set))
            {
              printf("Chamber pressure not set. Aborted.\n");
              break;
            }

            equil->properties.P = case_list[i].pressure;
                        
            print_propellant_composition(equil);
            TIME(if (equilibrium(equil, HP)) break,
                 CHAMBER_MSG);
            
            print_product_properties(equil, 1);
            print_product_composition(equil, 1);
            break;

        case FROZEN_PERFORMANCE:

            if (!(case_list[i].pressure_set))
            {
              printf("Chamber pressure not set. Aborted.\n");
              break;
            }
            else if (!(case_list[i].exit_condition_set))
            {
              printf("Exit condition not set. Aborted.\n");
              break;
            }
            
            copy_equilibrium(frozen, equil);
            
            equil->properties.T = case_list[i].temperature;
            equil->properties.P = case_list[i].pressure;
            
            print_propellant_composition(frozen);

            TIME(if (equilibrium(frozen, HP)) break,
                 CHAMBER_MSG);

            TIME(frozen_performance(frozen, case_list[i].exit_cond_type,
                                    case_list[i].exit_condition),
                 FROZEN_MSG);
            
            print_product_properties(frozen, 3);
            print_performance_information(frozen, 3);
            print_product_composition(frozen, 3);
            
          break;
        case EQUILIBRIUM_PERFORMANCE:

            
            if (!(case_list[i].pressure_set))
            {
              printf("Chamber pressure not set. Aborted.\n");
              break;
            }
            else if (!(case_list[i].exit_condition_set))
            {
              printf("Exit condition not set. Aborted.\n");
              break;
            }
            
            copy_equilibrium(shifting, equil);

            equil->properties.T = case_list[i].temperature;
            equil->properties.P = case_list[i].pressure;
            
            print_propellant_composition(shifting);
            
            TIME(if (equilibrium(shifting, HP)) break,
                 CHAMBER_MSG);

            TIME(equilibrium_performance(shifting, case_list[i].exit_cond_type,
                                         case_list[i].exit_condition),
                 EQUILIBRIUM_MSG);

            print_product_properties(shifting, 3);
            print_performance_information(shifting, 3);
            print_product_composition(shifting, 3);
            
            break;
      }
      i++;
    }
    free (equil);
    free (frozen);
    free (shifting);
    
  }
  
  free (propellant_list);
  free (thermo_list);

  if (errorfile != stderr)
    fclose (errorfile);

  if (outputfile != stdout)
    fclose (outputfile);
  
  return 0;

}
