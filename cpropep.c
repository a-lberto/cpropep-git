/* cpropep.c  -  Calculation of Complex Chemical Equilibrium           */
/* Copyright (C) 2000                                                  */
/* Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                      */
/* Contribution of                                                     */
/* Mark Pinese <ida.pinese@bushnet.qld.edu.au>                         */

/* This program is free software; you can redistribute it and/or modify*/
/* it under the terms of the GNU General Public License as published by*/
/* the Free Software Foundation; either version 2 of the License, or   */
/* (at your option) any later version.                                 */

/* This program is distributed in the hope that it will be useful,     */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of      */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */ 
/* GNU General Public License for more details.                        */

/* You should have received a copy of the GNU General Public License   */
/* along with this program; if not, write to the Free Software         */
/* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.           */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
/* #include <time.h> */

#include "print.h"
#include "equilibrium.h"
#include "load.h"
#include "libnum.h"

#include "performance.h"

#include "getopt.h"


#define version "1.0"
#define date "27/02/2000"


/* global variable containing the information about chemical species */
extern propellant_t  *propellant_list;
extern thermo_t	     *thermo_list;

extern double conc_tol;
extern double conv_tol;
extern int    iteration_max;


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

void info(char **argv) {
  printf("Try `%s -h' for more information\n", argv[0]);
}

void usage(void)
{
  printf("Option list\n");
  printf("-f file \t Input file\n");
  printf("-v num  \t Set the verbose level\n");
  printf("-p      \t Print the propellant list\n");
  printf("-t      \t print the combustion product list\n");

}


int load_input(FILE *fd, equilibrium_t *e, problem_t *p)
{ 
  double tmp1, tmp2;

  int sp;
  double m;
 
  int section = 0;

  char buffer[128];
  char tmp[10];
  char num[10];
  char qt[10];
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
            if ( strncmp(buffer, "TP", 2) == 0)
              *p = TP;
            else if (strncmp(buffer, "HP", 2) == 0)
              *p = HP;
            else
            {
              printf ("Unknown option.\n");
              break;
            }
            
            sscanf(buffer, "%s %lf %lf", tmp, &tmp1, &tmp2);
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
            //printf("Comments\n");
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
  int c;
  int v = 0;
  char filename[FILENAME_MAX];
  FILE *fd = NULL;
  equilibrium_t *equil;

  int thermo_loaded     = 0;
  int propellant_loaded = 0;

  problem_t p;
  
  /* allocate memory to hold data */
  if (mem_alloc())
    return 1;


  while (1)
  {
    c = getopt(argc, argv, "aphtf:v:");

    if (c == EOF)
      break;
    
    switch (c)
    {
      case 'a':
          printf("Got an a\n");
          break;

      case 'p':
          if (!propellant_loaded)
          {
            load_propellant ("propellant.dat");
            propellant_loaded = 1;
          }
          print_propellant_list();
          return 0;
          
      case 'h':
          usage();
          return 0;
          
      case 'f':
          if (strlen(optarg) > 128)
          {
            printf("Filename too long!\n");
            break;
          }
          strncpy(filename, optarg, 128);
          if ((fd = fopen(filename, "r")) == NULL )
            return 1;
          break;
          
      case 't':
          if (!thermo_loaded)
          {
            load_thermo ("thermo.dat");
            thermo_loaded = 1;
          }
          print_thermo_list();
          return 0;
          
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
          return 0;
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
  
  welcome_message();
  
  if (fd != NULL)
  {
    equil = (equilibrium_t *) malloc ( sizeof (equilibrium_t) );
    initialize_equilibrium(equil);
    load_input(fd, equil, &p);
    fclose(fd);
    set_verbose(equil, v);

    frozen_performance(equil, 1);
    //equilibrium(equil, p);

    dealloc_equillibrium (equil);
  }
 
  free (propellant_list);
  free (thermo_list);

  return 0;

}


