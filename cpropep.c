/* cpropep.c  -  Calculation of Complex Chemical Equilibrium           */
/* $Id: cpropep.c,v 1.20 2000/06/20 02:15:12 antoine Exp $ */
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
#define date    "23/05/2000"

#define CHAMBER_MSG     "Time spent for computing chamber equilibrium"
#define FROZEN_MSG      "Time spent for computing frozen performance"
#define EQUILIBRIUM_MSG "Time spent for computing equilibrium performance"

#define TIME(function, msg) timer = clock(); function;\
                            fprintf(outputfile,\
                                    "%s: %f s\n\n", msg,\
                                    (double)(clock() - timer)/CLOCKS_PER_SEC)

#define THERMO_FILE "new_thermo.dat"
#define PROPELLANT_FILE "propellant.dat"


//#undef TIME
//#define TIME(function, msg) function;

#define MAX_CASE 10


typedef enum _p
{
  SIMPLE_EQUILIBRIUM,
  FIND_FLAME_TEMPERATURE,
  MULTIPLE_FLAME_TEMPERATURE,
  FROZEN_PERFORMANCE,
  EQUILIBRIUM_PERFORMANCE,
  ALL_PERFORMANCE,
  MULTIPLE_PERFORMANCE
} p_type;

char case_name[][80] = {
  "Fixed pressure-temperature equilibrium",
  "Fixed enthalpy-pressure equilibrium - adiabatic flame temperature",
  "Multiple adiabatic flame temperature problem",
  "Frozen equilibrium performance evaluation",
  "Shifting equilibrium performance evaluation",
  "Frozen and shifting equilibrium performance evaluation",
  "Multiple performance evaluation problem"
};

typedef struct _case_t
{
  p_type p;
  double arg[6];
} case_t;


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
  printf("-q num  \t Print information about propellant number num\n");
  printf("-t      \t Print the combustion product list\n");
  printf("-u num  \t Print information about product number num\n");
  printf("-h      \t Print help\n");
  printf("-i      \t Print informations\n");
}


//int load_input(FILE *fd, equilibrium_t *e, p_type *p, double *pe)
int load_input(FILE *fd, equilibrium_t *e, case_t *t, double *pe)
{ 
  double m;//tmp1, tmp2;
  
  int sp;
  int section = 0;

  int n_case = 0;
  
  char buffer[128], tmp[10], num[10], qt[10];
  
  char *bufptr;

  char unit; /* unit of measure of the quantity */

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
            section = 1;
          else
          { 
            if (strncmp(buffer, "TP", 2) == 0)
              t[n_case].p = SIMPLE_EQUILIBRIUM;
            else if (strncmp(buffer, "HP", 2) == 0)
              t[n_case].p = FIND_FLAME_TEMPERATURE;
            else if (strncmp(buffer, "MHP", 2) == 0)
              t[n_case].p = MULTIPLE_FLAME_TEMPERATURE;
            else if (strncmp(buffer, "FR", 2) == 0)
              t[n_case].p = FROZEN_PERFORMANCE;
            else if (strncmp(buffer, "EQ", 2) == 0)
              t[n_case].p = EQUILIBRIUM_PERFORMANCE;
            else if (strncmp(buffer, "PE", 2) == 0)
              t[n_case].p = ALL_PERFORMANCE;
            else if (strncmp(buffer, "MP", 2) == 0)
              t[n_case].p = MULTIPLE_PERFORMANCE;
            else
            {
              printf ("Unknown option.\n");
              break;
            }

            sscanf(buffer, "%s %lf %lf %lf %lf %lf %lf",
                   tmp,
                   t[n_case].arg,
                   t[n_case].arg+1,
                   t[n_case].arg+2,
                   t[n_case].arg+3,
                   t[n_case].arg+4,
                   t[n_case].arg+5);

            *pe = t[n_case].arg[2];
            //sscanf(buffer, "%s %lf %lf %lf", tmp, &tmp1, &tmp2, pe);
            set_state(e, t[n_case].arg[0], t[n_case].arg[1]);
            n_case++;
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
  int i, c, v = 0;
  char filename[FILENAME_MAX];
  FILE *fd = NULL;

  equilibrium_t *equil, *frozen, *shifting; 
  
  int thermo_loaded     = 0;
  int propellant_loaded = 0;

  double exit_pressure;

  //int temp;

  clock_t timer;

  int param;
  
//  p_type p;

  case_t case_list[MAX_CASE];
  for (i = 0; i < MAX_CASE; i++)
    case_list[i].p = -1;

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
    
    fclose(fd);
    global_verbose = v;
    //set_verbose(equil, v);

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
            set_state(equil, case_list[i].arg[0], case_list[i].arg[1]);
            
            print_propellant_composition(equil);
            TIME(if (equilibrium(equil, TP)) break,
                 CHAMBER_MSG);

            print_product_properties(equil, 1);
            print_product_composition(equil, 1);
            break;
        case FIND_FLAME_TEMPERATURE:
            set_state(equil, case_list[i].arg[0], case_list[i].arg[1]);
            
            print_propellant_composition(equil);
            TIME(if (equilibrium(equil, HP)) break,
                 CHAMBER_MSG);
            
            print_product_properties(equil, 1);
            print_product_composition(equil, 1);
            break;
        case MULTIPLE_FLAME_TEMPERATURE:
            print_propellant_composition(equil);

            for (param = case_list[i].arg[0];
                 param < case_list[i].arg[2];
                 param += case_list[i].arg[1])
            {
              set_state(equil, 0, param);

              equilibrium(equil, HP);
              
              print_product_properties(equil, 1);
              print_product_composition(equil, 1);
            }
            break;
              
        case FROZEN_PERFORMANCE:

            copy_equilibrium(frozen, equil);
            
            exit_pressure = case_list[i].arg[2];
            set_state(frozen, case_list[i].arg[0], case_list[i].arg[1]);

            print_propellant_composition(frozen);

            TIME(if (equilibrium(frozen, HP)) break,
                 CHAMBER_MSG);

            TIME(frozen_performance(frozen, exit_pressure),
                 FROZEN_MSG);
            
            print_product_properties(frozen, 3);
            print_performance_information(frozen, 3);
            print_product_composition(frozen, 3);
            
          break;
        case EQUILIBRIUM_PERFORMANCE:
            copy_equilibrium(shifting, equil);
            
            exit_pressure = case_list[i].arg[2];
            set_state(shifting, case_list[i].arg[0], case_list[i].arg[1]);
            
            print_propellant_composition(shifting);
            
            TIME(if (equilibrium(shifting, HP)) break,
                 CHAMBER_MSG);

            TIME(equilibrium_performance(shifting, exit_pressure),
                 EQUILIBRIUM_MSG);

            print_product_properties(shifting, 3);
            print_performance_information(shifting, 3);
            print_product_composition(shifting, 3);
            
            break;
        case ALL_PERFORMANCE:
/*
            copy_equilibrium(frozen, equil);
            copy_equilibrium(shifting, equil);
            
            exit_pressure = case_list[i].arg[2];
            set_state(frozen, case_list[i].arg[0], case_list[i].arg[1]);
            set_state(shifting, case_list[i].arg[0], case_list[i].arg[1]);
            
            print_propellant_composition(equil);

            TIME(if (equilibrium(equil, HP)) break,
               CHAMBER_MSG);
            printf("--- Chamber equilibrium properties ---\n");
            derivative(equil, &deriv);
            print_product_properties(equil);
            print_derivative_results(&deriv);
            print_product_composition(equil);
            TIME(frozen_performance(equil, &performance, exit_pressure),
                 FROZEN_MSG);
            TIME(equilibrium_performance(equil, exit_equil, &performance,
                                         exit_pressure),
                 EQUILIBRIUM_MSG);
            print_performance_information(&performance);
            
            printf("--- Exit equilibrium properties ---\n");
            print_product_properties(exit_equil);
            print_product_composition(exit_equil);
*/
            break;
        case MULTIPLE_PERFORMANCE:
/*
            set_state(equil, case_list[i].arg[0], case_list[i].arg[1]);
            
            print_propellant_composition(equil);
            TIME(if (equilibrium(equil, HP)) break,
               CHAMBER_MSG);
            printf("--- Chamber equilibrium properties ---\n");
            derivative(equil, &deriv);
            print_product_properties(equil);
            print_derivative_results(&deriv);
            print_product_composition(equil);
            for (exit_pressure = case_list[i].arg[2];
                 exit_pressure < case_list[i].arg[4];
                 exit_pressure += case_list[i].arg[3])
            {
              TIME(frozen_performance(equil, &performance, exit_pressure),
                   FROZEN_MSG);
              TIME(equilibrium_performance(equil, exit_equil, &performance,
                                           exit_pressure),
                   EQUILIBRIUM_MSG);
              print_performance_information(&performance);
            }
            break;
*/
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
