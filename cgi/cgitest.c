/*
    cgitest.c - Testprogram for cgi.o
    Copyright (c) 1998,9 by Martin Schulze <joey@infodrom.north.de>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 */

/*
 * Compile with: cc -o cgitest cgitest.c -lcgi
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cgi.h>


#include "equilibrium.h"
#include "performance.h"
#include "thermo.h"

#include "load.h"
#include "print.h"


s_cgi *cgi;

#ifdef SOURCEFORGE
#define THERMO_PATH "thermo.dat"
#define PROPELLANT_PATH "propellant.dat"
#else
#define THERMO_PATH "/home/antoine/projets/rocket/rocketworkbench/cpropep/thermo.dat"
#define PROPELLANT_PATH "/home/antoine/projets/rocket/rocketworkbench/cpropep/propellant.dat"
#endif

void init_equil(void)
{
  load_thermo (THERMO_PATH);
  load_propellant (PROPELLANT_PATH); 
}

void destroy(void)
{
  free (propellant_list);
  free (thermo_list);
}

int eval_cgi(equilibrium_t *e, double *exit_condition,
             exit_condition_t *exit_type)
{
  int    mol = 0;
  char   buffer[64];
  double tmp;
  double tmp1;
  
  printf ("<h1>Results</h1>\n\n");
    
  strncpy(buffer, cgiGetValue(cgi, "temp"), 64);
  
  if ( ((tmp = atof(buffer)) != 0) && tmp > 298.15 && tmp < 6000 )
    e->properties.T = tmp;
  else 
    return -1;
  
  strncpy(buffer, cgiGetValue(cgi, "pressure"), 64);

  if ((tmp = atof(buffer)) != 0)
    e->properties.P = tmp;
  else
    return -1;

  strncpy(buffer, cgiGetValue(cgi, "exit_cond_type"), 64);

  if (!(strcmp(buffer, "Pressure")))
    *exit_type = PRESSURE;
  else if (!(strcmp(buffer, "Supersonic aera ratio")))
    *exit_type = SUPERSONIC_AREA_RATIO;
  else if (!(strcmp(buffer, "Subsonic aera ratio")))
    *exit_type = SUBSONIC_AREA_RATIO;
  else
    return -1;

  strncpy(buffer, cgiGetValue(cgi, "exit_cond"), 64);

  if ((tmp = atof(buffer)) != 0)
    *exit_condition = tmp;
  else
    return -1;

    
  if (strncmp("mol", cgiGetValue(cgi, "select"), 3) == 0)
    mol = 1;
  
  strncpy(buffer, cgiGetValue(cgi, "qt1"), 32);
  tmp = atof(buffer);
  
  strncpy(buffer,  cgiGetValue(cgi, "ml1"), 32);
  tmp1 = atof(buffer);
  
  if ((tmp != 0) && (tmp1 != 0))
  {
    if (mol)
	    add_in_propellant(e, tmp1, tmp);
    else
	    add_in_propellant(e, tmp1, GRAM_TO_MOL(tmp, tmp1));
  }
  
  strncpy(buffer, cgiGetValue(cgi, "qt2"), 32);
  tmp = atof(buffer);

  strncpy(buffer, cgiGetValue(cgi, "ml2"), 32);
  tmp1 = atof(buffer);
  
  if ((tmp != 0) && (tmp1 != 0))
  {
    if (mol)
	    add_in_propellant(e, tmp1, tmp);
    else
	    add_in_propellant(e, tmp1, GRAM_TO_MOL(tmp, tmp1));
  }
  
  strncpy(buffer, cgiGetValue(cgi, "qt3"), 32);
  tmp = atof(buffer);
  
  strncpy(buffer, cgiGetValue(cgi, "ml3"), 32);
  tmp1 = atof(buffer);
  
  if ((tmp != 0) && (tmp1 != 0))
  {
    if (mol)
	    add_in_propellant(e, tmp1, tmp);
    else
	    add_in_propellant(e, tmp1, GRAM_TO_MOL(tmp, tmp1));
  }

  compute_density(&(e->propellant));
  
  return 0;
}

int main (int argc, char **argv, char **env)
{
  short i;
  
  char *path_info = NULL;  
 
  equilibrium_t *equil, *frozen, *shifting; 
  
  char buffer[32];

  double value;
  double exit_cond;
  exit_condition_t exit_cond_type;
  
  int val;
  problem_t P = TP;

  errorfile = stderr;
  outputfile = stdout;
    
  cgiDebug(0, 0);
  cgi = cgiInit();
  
  init_equil();

  global_verbose = 0;
  
  path_info = getenv("PATH_INFO");

  if (path_info) 
  {
    if (!strcmp(path_info, "/list")) 
    {
      cgiHeader();
	    printf("<pre>");
	    print_propellant_list();
	    printf("</pre>");
	    
    }
    else if (!strcmp(path_info, "/list_thermo"))
    {
      cgiHeader();
	    printf("<pre>");
	    print_thermo_list();
	    printf("</pre>");
      
    }
    else if (!strcmp(path_info, "/search_thermo"))
    {
      cgiHeader();
      printf("<pre>");
      thermo_search( cgiGetValue(cgi, "formula"));
      printf("</pre>");
      
      
    }
    else if (!strcmp(path_info, "/search_prop"))
    {
      cgiHeader();
      printf("<pre>");
      propellant_search( cgiGetValue(cgi, "name"));
      //thermo_search( cgiGetValue(cgi, "formula"));
      printf("</pre>");

    } 
    else if (!strcmp(path_info, "/equil")) 
    {
      cgiHeader();

      equil = (equilibrium_t *) malloc (sizeof (equilibrium_t));
      initialize_equilibrium(equil);
      
      frozen   = (equilibrium_t *) malloc (sizeof(equilibrium_t)*3);
      shifting = (equilibrium_t *) malloc (sizeof(equilibrium_t)*3);
      for (i = 0; i < 3; i++)
      {
        initialize_equilibrium(frozen + i);
        initialize_equilibrium(shifting + i);
      }
      
	    if (eval_cgi(equil, &exit_cond, &exit_cond_type))
        printf("<b>Error</b>");
	    else
	    {
        list_element(equil);
        list_product(equil);
      
        printf("<pre>");
        

        if (strncmp("Find", cgiGetValue(cgi, "type"), 4) == 0)
        {
          print_propellant_composition(equil);
          P = HP;
          equilibrium(equil, P);
          print_product_properties(equil, 1);
          print_product_composition(equil, 1);
          
        }
        else if (strncmp("Froz", cgiGetValue(cgi, "type"), 4) == 0)
        {
          copy_equilibrium(frozen, equil);
          
          print_propellant_composition(frozen);
          equilibrium(frozen, HP);
            
          frozen_performance(frozen, exit_cond_type, exit_cond);

          print_product_properties(frozen, 3);
          print_performance_information(frozen, 3);
          print_product_composition(frozen, 3);
          
        }
        else if (strncmp("Shifting", cgiGetValue(cgi, "type"), 8) == 0)
        {
          copy_equilibrium(shifting, equil);
          
          print_propellant_composition(shifting);

          equilibrium(shifting, HP);
          shifting_performance(shifting, exit_cond_type, exit_cond);

          print_product_properties(shifting, 3);
          print_performance_information(shifting, 3);
          print_product_composition(shifting, 3);
        }
        else
        {
          print_propellant_composition(equil);
          equilibrium(equil, P);
          print_product_properties(equil, 1);
          print_product_composition(equil, 1);
          
        }
          
        printf("</pre>");
	    }
      
      free (equil);
      free (frozen);
      free (shifting);
      
    }
    else if (!strcmp(path_info, "/prop"))
    {
      cgiHeader();
      printf("<pre>");
      strncpy(buffer, cgiGetValue(cgi, "propellant"), 32);

      if ( print_propellant_info (atoi (buffer)))
        printf("Request out of range\n");
     
      printf("</pre>");
      
    }
    else if (!strcmp(path_info, "/prod"))
    {
      
      cgiHeader();
      printf("<pre>");
      strncpy(buffer, cgiGetValue(cgi, "product"), 32);
      val = atoi (buffer);
      
      if ( print_thermo_info (val))
        printf("Request out of range\n");
      else
      {
        strncpy(buffer, cgiGetValue(cgi, "temp"), 32);
        value = atof(buffer);

        printf("Thermodynamics properties at temperature %.2f K\n", value);
        printf("Enthalpy:      % f J/mol\n", enthalpy_0(val, value)*R*value);
        printf("Entropy:       % f J/(mol K)\n", entropy_0(val, value)*R);
        printf("Specific heat: % f J/(mol K)\n",
               specific_heat_0(val, value)*R);
      
      }
      printf("</pre>");

    }
    else 
    {
	    cgiHeader();
	    printf("<br>Bad queries<br>");
    }
  } 
  else 
  {
    cgiHeader();
    printf("<br>Bad PATH_INFO<br>");
  }
  
  printf ("\n<hr>\n</body>\n</html>\n");

  
  destroy();
  return 0;
}

