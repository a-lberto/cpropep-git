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

#include "performance.h"
#include "equilibrium.h"
#include "load.h"
#include "print.h"

extern propellant_t  *propellant_list;
extern thermo_t	     *thermo_list;

extern FILE * errorfile;
extern FILE * outputfile;

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

int eval_cgi(equilibrium_t *e)
{
  int    mol = 0;
  char   buffer[32];
  double tmp;
  double tmp1;
  
  printf ("<h1>Results</h1>\n\n");
    
  strncpy(buffer, cgiGetValue(cgi, "temp"), 32);
  
  if ( ((tmp = atof(buffer)) != 0) && tmp > 298.15 && tmp < 6000 )
    e->T = tmp;
  else 
    return -1;
  
  strncpy(buffer, cgiGetValue(cgi, "pressure"), 32);

  if ((tmp = atof(buffer)) != 0)
    e->P = tmp;
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
  return 0;
}

int main (int argc, char **argv, char **env)
{
  char *path_info = NULL;  
  equilibrium_t *equil, *exit_equil;
  performance_t performance;
  
  char buffer[32];

  double value;
  int val;
  problem_t P = TP;

  errorfile = stderr;
  outputfile = stdout;
  
  performance.frozen_ok = false;
  performance.equilibrium_ok = false;
  
  cgiDebug(0, 0);
  cgi = cgiInit();
  
  init_equil();
  
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
      
	    equil = (equilibrium_t *) malloc ( sizeof (equilibrium_t) );
	    initialize_equilibrium(equil);
      exit_equil = (equilibrium_t *) malloc ( sizeof (equilibrium_t) );
	    initialize_equilibrium(exit_equil);
      
      
	    if (eval_cgi(equil))
        printf("<b>Error</b>");
	    else
	    {
        printf("<pre>");
        set_verbose(equil, 0);
        set_verbose(exit_equil, 0);

        if (strncmp("Find", cgiGetValue(cgi, "type"), 4) == 0)
        {
          print_propellant_composition(equil);
          P = HP;
          equilibrium(equil, P);
          print_product_properties(equil);
          print_product_composition(equil);
          
        }
        else if (strncmp("Froz", cgiGetValue(cgi, "type"), 4) == 0)
        {
          print_propellant_composition(equil);
          equilibrium(equil, HP);

          printf("--- Chamber equilibrium properties ---\n");
          print_product_properties(equil);
          print_product_composition(equil);
            
          frozen_performance(equil, &performance, 1);
          print_performance_information(&performance);
        }
        else if (strncmp("Shifting", cgiGetValue(cgi, "type"), 8) == 0)
        {
          print_propellant_composition(equil);
          equilibrium(equil, HP);

          printf("--- Chamber equilibrium properties ---\n");
          print_product_properties(equil);
          print_product_composition(equil);
          
          equilibrium_performance(equil, exit_equil, &performance, 1);
          print_performance_information(&performance);

           printf("--- Exit equilibrium properties ---\n");
           print_product_properties(exit_equil);
           print_product_composition(exit_equil);
        }
        else
        {
          print_propellant_composition(equil);
          equilibrium(equil, P);
          print_product_properties(equil);
          print_product_composition(equil);
          
        }
          
        printf("</pre>");
	    }
      
	    dealloc_equilibrium (equil);
      dealloc_equilibrium (exit_equil);
      free (equil);
      free (exit_equil);
      
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

