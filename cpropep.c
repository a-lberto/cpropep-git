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

#include "equilibrium.h"
#include "load.h"
#include "libnum.h"

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

int main(int argc, char *argv[])
{
  
  int i;

  equilibrium_t *equil;

  //for (i = 0; i < argc; i++)
  //{
  //  if ( argv[i][0] == '-' )
  //  {
  //    switch (argv[i][1])
  //    {
  //    case 'a':
  //	printf("%s\n", (argv[i] + 1) );
  //	break;
  //    case 'c':
  //	printf("%s\n", (argv[i] + 1) );
  //	break;
  //    default:
  //	printf("Unknown option\n");
  //   }
  //  }
  //}
  //printf("%s\n",argv[i]);
    
  equil = (equilibrium_t *) malloc ( sizeof (equilibrium_t) );
  initialize_equilibrium(equil);

  /* set the state for equil, temperature, pressure */
  set_state(equil, 2598, 6);

  /* allocate memory to hold data */
  if (mem_alloc())
    return 1;
   

  initialize_product(equil->p);

  load_thermo ("thermo.dat");
  load_propellant ("propellant.dat");


  welcome_message();

  //add_in_propellant(equil, 685, GRAM_TO_MOL(79.36, 685) ); // O2
  //add_in_propellant(equil, 458, GRAM_TO_MOL(10, 458) ); // H2

  //add_in_propellant(equil, 963, 1); // UDMH
  //add_in_propellant(equil, 436, 1); // HYDRAZINE N2H4
  //add_in_propellant(equil, 378, 2); // F2

  add_in_propellant(equil, 766, GRAM_TO_MOL(70, 766) ); // KClO4
  //add_in_propellant(equil, 788, 2); // HTPB

  //add_in_propellant(equil, 765, 0.742); // KNO3
  add_in_propellant(equil, 840, GRAM_TO_MOL(30, 840)); // table sugar H22C12O11
  //add_in_propellant(equil, 210, 1.25); // C

  set_verbose(equil, 1);
  
  equilibrium(equil);


  //print_thermo_list();
  //print_propellant_list();
    


  dealloc_equillibrium (equil);
 
  free (propellant_list);
  free (thermo_list);

  return 0;
  
}



