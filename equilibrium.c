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


/* global variable containing the information about chemical species */
propellant_t	*propellant_list;
thermo_t	*thermo_list;


double conc_tol = 1e-8;
double conv_tol = 0.5e-5;
int iteration_max = 35;

/****************************************************************
VARIABLE: Contain the molar mass of element by atomic number
          molar_mass[0] contain hydrogen and so on.
*****************************************************************/
const float molar_mass[100] = { 
  1.00794, 4.00260, 6.941, 9.01218, 10.811, 12.011,
  14.0067, 15.9994, 18.99840, 20.11797, 22.99977, 24.305, 
  26.98154, 28.0855, 30.97376, 32.066, 35.4527, 39.948, 
  39.0983, 40.078, 44.9559, 47.88, 50.9415, 51.996, 54.938, 
  55.847, 58.9332, 58.6934, 63.546, 65.39, 69.723, 72.61, 
  74.9216,
  // we should correct the value
  78.96,79.916, 83.80, 85.48, 87.63, 88.91, 91.22, 92.91, 95.95,
  99.,101.1, 102.91, 106.4, 107.88, 112.41, 114.82, 118.7, 121.76,
  127.61, 126.91, 131.3, 132.91, 137.36, 138.92, 140.13, 140.91,
  144.27,147., 150.35, 152., 157.26, 158.93, 162.51,164.94,167.27,
  168.94,173.04, 174.99, 178.50, 180.95, 183.86, 186.22, 190.2,
  192.2,195.09, 197., 220.61, 204.39, 207.21, 208.99, 210., 210.,
  222.,2.014,226.,92.906,232.,231.,238.,237.,237.,12.011,9.013,
  10.82,24.32,26.98, 253.0 };


/****************************************************************
VARIABLE: Contain the symbol of the element in the same way as
          for the molar mass.

COMMENTS: It is use in the loading of the data file to recognize
          the chemical formula.
*****************************************************************/
const char symb[N_SYMB][3] = {
  "H ","HE","LI","BE","B ","C ","N ","O ",
  "F ","NE","NA","MG","AL","SI","P ","S ","CL","AR","K ","CA",
  "SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN","GA","GE",
  "AS","SE","BR","KR","RB","SR","Y ","ZR","NB","MO","TC","RU",
  "RH","PD","AG","CD","IN","SN","SB","TE","I ","XE","CS","BA",
  "LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER",
  "TM","YB","LU","HF","TA","W ","RE","OS","IR","PT","AU","HG","TL",
  "PB","BI","PO","AT","RN","FR","RA","AC","TH","PA","U ","NP",
  "U6","U5","U1","U2","U3","U4","FM",
  "E ", "D " }; /* the E stand for electron and D for deuterium*/


/***************************************************************
MACRO: molar gaz constant in J/(mol K)
****************************************************************/
const float  R = 8.3143; 


double enthalpy(int sp, float T)
{
  thermo_t s = *(thermo_list + sp);
  double val;
  int i;
     
  for (i = 0; i < 4; i++)
  {
    if ((T >= s.range[i][0]) && (T < s.range[i][1])) 
    {
      /* parametric equation for dimentionless enthalpy */
      val = -s.param[i][0]*pow(T, -2) + s.param[i][1]*pow(T, -1)*log(T) + s.param[i][2] + s.param[i][3]*T/2 + s.param[i][4]*pow(T, 2)/3 + s.param[i][5]*pow(T, 3)/4 + s.param[i][6]*pow(T, 4)/5 + s.param[i][7]/T;

      return val; /* dimensionless enthalpy */
    }
  }
  
  printf("Error: temperature out of range for %d: %s\n", sp, s.name);
  return 0;
}


double entropy(int sp, float T)
{
  thermo_t s = thermo_list[sp];
  double val;
  int i;
  
  for (i = 0; i < 4; i++)
  {
    if ((T >= s.range[i][0]) && (T < s.range[i][1]))
    {	
      /* parametric equation for dimentionless entropy */
      val = -s.param[i][0]*pow(T, -2)/2 - s.param[i][1]*pow(T, -1)
	+ s.param[i][2]*log(T) + s.param[i][3]*T +s.param[i][4]*pow(T, 2)/2
	+ s.param[i][5]*pow(T, 3)/6 + s.param[i][6]*pow(T, 4)/4 
	+ s.param[i][8];
      
      return val; /* dimensionless entropy */
    }
  }
  
  printf("Error: temperature out of range for %d: %s\n", sp, s.name);
  return 0;
}

double specific_heat(int sp, float T)
{
  thermo_t s = thermo_list[sp];
  double val;
  int i;
     
  for (i = 0; i < 4; i++)
  {
    if ((T >= s.range[i][0]) && (T < s.range[i][1])) 
    {	
      /* parametric equation for dimentionless entropy */
      val = s.param[i][0]*pow(T, -2) + s.param[i][1]*pow(T, -1)
	+ s.param[i][2] + s.param[i][3]*T +s.param[i][4]*pow(T, 2)
	+ s.param[i][5]*pow(T, 3) + s.param[i][6]*pow(T, 4);
      
      return val; /* dimensionless */
    }
  }
  
  printf("Error: temperature out of range for %d: %s\n", sp, s.name);
  return 0;
}

/* 0 if out of range, 1 if ok */
int temperature_check(int sp, float T)
{
  int i;
  thermo_t s = thermo_list[sp];
  for (i = 0; i < 4; i++)
  {
    if ((T >= s.range[i][0]) && (T < s.range[i][1]))
      return 1;
  }
    return 0;
}

/* unuse function */
double delta_enthalpy(int sp, float T)
{
  /* delta henthalpy in J/mol */
  return enthalpy(sp, T) - (thermo_list + sp)->heat;
}


double gibbs0(int sp, float T)
{
  return enthalpy(sp, T) - entropy(sp, T);
}

/* J/mol T is in K, P is in atm */
/* uo = HO(T) - S(T)  (from the nasa book) */
double gibbs(int sp, state_t st, double nj, double n, float T, float P)
{

  double g;

  if ( (st == GAS) && (nj == 0.0))
    return 0;

  switch (st)
  {
  case GAS:    
    g = gibbs0(sp, T) + log(nj/n) + log(P);
    break;
  case CONDENSED:
    g = gibbs0(sp, T);
    break;
  default:
    g = 0;
  }  
  return g;
}


double propellant_molar_mass(int molecule)
{     
  int i;
  double ans = 0;
  for (i = 0; i < 6; i++)
    ans += (propellant_list + molecule)->coef[i] * molar_mass[(propellant_list + molecule)->elem[i]];
  return ans;
}

/* kJ/mol */
double heat_of_formation(int molecule)
{
  return (propellant_list + molecule)->heat * 4.1868 * 
    propellant_molar_mass(molecule) / 1000;
}



int print_thermo_info(int sp)
{
  int i;
  int j;
  
  thermo_t s = *(thermo_list + sp);
  
  printf("---------------------------------------------\n");
  printf("Name: \t\t\t%s\n", s.name);
  printf("Comments: \t\t%s\n", s.comments);
  printf("Id: \t\t\t%s\n", s.id);
  printf("Chemical formula:\t");
  
  for (i = 0; i < 5; i++)
  {
    if (!(s.coef[i] == 0))
      printf("%d%s", s.coef[i], symb[ s.elem[i]]);
  }
  printf("\n");
  printf("State:\t\t\t");
  switch (s.state)
  {
  case GAS:
    printf("GAZ\n");
    break;
  case CONDENSED:
    printf("CONDENSED\n");
    break;
  default:
    printf("UNRECOGNIZE\n");
  }
  
  printf("\n");
  printf("Molecular weight: \t\t%f g/mol\n", s.weight);
  printf("Heat of formation at 298.15 K : %f J/mol\n", s.heat);
  printf("HO(298.15) - HO(0): \t\t%f J/mol\n", s.dho);
  printf("Number of temperature range: %d\n\n", s.nint);
  
  for (i = 0; i < s.nint; i++)
  {
    printf("Interval: %f - %f \n\n", s.range[i][0], s.range[i][1]);
    for (j = 0; j < 9; j++)
      printf("%Le ", s.param[i][j]);
    printf("\n");
  }
  printf("---------------------------------------------\n");
  return 0;
}

int thermo_search(const char *str)
{
  int i;    

  /* Note - must be changed in the future */

  for (i = 0; i < MAX_THERMO; i++)
  {
    if (!(strncmp(str, (thermo_list + i)->name, strlen(str))))
      return i;
  }
  return -1;
}


int print_thermo_list(void)
{
  int i;
  for (i = 0; i < MAX_THERMO; i++)
    printf("%d   %s\n", i, (thermo_list + i)->name);
 
  return 0;
}

int print_propellant_list(void)
{
  int i;
  for (i = 0; i < MAX_PROPELLANT; i++)
    printf("%d   %s\n", i, (propellant_list + i)->name);
 
  return 0;
}

int print_condensed(product_t p)
{
  int i;
  for (i = 0; i < p.n[CONDENSED]; i ++)
    printf("%s ",  (thermo_list + p.species[CONDENSED][i])->name );
  printf("\n");
  return 0;
}

int print_gazeous(product_t p)
{
  int i;
  for (i = 0; i < p.n[GAS]; i++)
    printf("%s ", (thermo_list + p.species[GAS][i])->name );
  printf("\n");
  return 0;
}

int print_product_composition(equilibrium_t *e)
{
  int i;

  printf("%.4e mol of gazes\n", e->n);
  printf("\t  molar fraction  mol\n");
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    if (!(e->p->coef[GAS][i] == 0.0))
      printf("%s \t % .4e \t % .4e\n", 
	     (thermo_list + e->p->species[GAS][i])->name,
	     e->p->coef[GAS][i]/e->n,
	     e->p->coef[GAS][i]);
  }
  if (e->p->n[CONDENSED] > 0)
    printf("Condensed species (mol)\n");
  for (i = 0; i < e->p->n[CONDENSED]; i++)
  {
    printf("%s  % .4e\n", (thermo_list + e->p->species[CONDENSED][i])->name,
	   e->p->coef[CONDENSED][i]);
  }

  printf("Molar mass of product: %f g/mol\n", product_molar_mass(e));

  return 0;
}

int print_propellant_composition(equilibrium_t *e)
{
  int i, j;

  printf("Propellant composition\n");
  printf("Code %-35s mol    Mass (g)  Composition\n", "Name");
  for (i = 0; i < e->c->ncomp; i++)
  {
    printf("%d  %-35s %.4f %.4f ", e->c->molecule[i],
	   (propellant_list + e->c->molecule[i])->name,
	   e->c->coef[i], 
	   e->c->coef[i]*propellant_molar_mass( e->c->molecule[i] ) );

    printf("  ");
    /* print the composition */
    for (j = 0; j < 6; j++)
    {
      if (!((propellant_list + e->c->molecule[i])->coef[j] == 0))
	printf("%d%s ", (propellant_list + e->c->molecule[i])->coef[j],
	       symb[ (propellant_list + e->c->molecule[i])->elem[j] ]);
    }
    printf("\n");
  }
  printf("\n");
  return 0;

}


double product_molar_mass(equilibrium_t *e)
{
  int i;
  double mol = 0.0;
  double mass = 0.0;

  for (i = 0; i < e->p->n[GAS]; i++)
  {
    mol += e->p->coef[GAS][i]; 
    mass += e->p->coef[GAS][i]*(thermo_list + e->p->species[GAS][i])->weight;
    
  }

  for (i = 0; i < e->p->n[CONDENSED]; i++)
    mass += e->p->coef[CONDENSED][i]*(thermo_list + 
				      e->p->species[CONDENSED][i])->weight;

  return (mass/mol);
}

int list_element(equilibrium_t *e)
{
  int n = 0;
  int t = 0;
  int i, j, k;
  
  for (i = 0; i < e->c->ncomp; i++)
  {
    for (j = 0; j < 6; j++)
    {	       
      if (!(propellant_list[ e->c->molecule[i] ].coef[j] == 0))
      {
	t = propellant_list[ e->c->molecule[i] ].elem[j];
	for (k = 0; k <= n; k++)
	{
	  if (e->element[k] == t)
	    break;
	  if (k == n)
	  {
	    e->element[n] = t;
	    n++;
	    break;
	  }
	}
      }
    }
  }
  e->n_element = n;

  printf("%d different elements in the propellant\n", n);
  /* Print those elements */
  for (i = 0; i < n; i++)
    printf("%s ", symb[e->element[i]] );
  printf("\n");

  return n;
}

int list_product(equilibrium_t *e)
{
  int i, j, k;

  int n = 0;   // global counter (number of species found)
  int st;      // temporary variable to hold the state of one specie
  int ok = 1;
  
  for (j = 0; j < MAX_THERMO; j++)
  {
    for (k = 0; k < 5; k++)
    {
      if (!(thermo_list[j].coef[k] == 0))
      {
	for (i = 0; i < e->n_element; i++)
	{
	  if (e->element[i] == thermo_list[j].elem[k])
	    break;
	  else if (i == (e->n_element - 1) )
	    ok = 0;
	}    
	if (!ok)
	  break;
      }
    }
    if (ok) /* add to the list */
    {
      st = thermo_list[j].state;

      if (temperature_check(j, e->T) )  // if the molecule could 
	                             // exist at that temperature
      {
	e->p->species[st][ e->p->n[st] ] = j;
	e->p->n[st]++;
	n++;
      }

    }
    ok = 1;
  }

  // reallocate the momory to the minimum size
  for (i = 0; i < STATE_LAST; i++)
  {
    if (( e->p->species[i] = (int *) realloc( e->p->species[i], 
					   sizeof(int) * e->p->n[i])) == NULL)
      printf("Reallocation of memory failed\n");
    
    if (( e->p->coef[i] = (double *) realloc ( e->p->coef[i], 
					sizeof(double) * e->p->n[i])) == NULL)
      printf("Reallocation of memory failed\n");
  }

  /* initialize tho mol number to 0.1mol/(nb of gazeous species) */
  for (i = 0; i < e->p->n[GAS]; i++)
    e->p->coef[GAS][i] = 0.1/e->p->n[GAS];
 
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 0.1, equil->n ....


  /* initialize condensed to zero */
  for (i = 0; i < e->p->n[CONDENSED]; i++)
    e->p->coef[CONDENSED][i] = 0;

  
  printf("%d possible combustion product\n", n);
  printf("%d gazeous species\n", e->p->n[GAS]);
  print_gazeous(*(e->p));
  printf("%d condensed species\n", e->p->n[CONDENSED]);
  print_condensed(*(e->p));


  return n;
}



int mem_alloc(void)
{
  if ((thermo_list = 
       (thermo_t *) malloc(sizeof(thermo_t) * MAX_THERMO)) == NULL)
  {
    printf ("\n\nMemory allocation error with thermo_t thermo_list[%i], %d bytes required", 
	    MAX_THERMO, sizeof(thermo_t) * MAX_THERMO);
    return 1;
  }
  else if ((propellant_list = 
	    (propellant_t *) malloc(sizeof(propellant_t) * MAX_PROPELLANT))
	   == NULL)
  {
    printf ("\n\nMemory allocation error with propellant_t propellant_list[%i], %d bytes required", 
	    MAX_PROPELLANT, sizeof(propellant_t) * MAX_PROPELLANT);
    /* Clean up */
    free (thermo_list);
    return 1;
  }
  return 0;
}

int initialize_product(product_t *p)
{
  int i, j;

  for (i = 0; i < STATE_LAST; i++)
    p->n[i] = 0;

  // allocate the memory
  for (i = 0; i < STATE_LAST; i++)
    p->species[i] = (int *) malloc (sizeof(int) * 300);

  for (i = 0; i < STATE_LAST; i++)
    p->coef[i] = (double *) malloc (sizeof(double) * 300);

  p->isalloc = 1;

  // initialize the list to -1
  for (i = 0; i < 300; i++)
  {
    for (j = 0; j < STATE_LAST; j++)
      p->species[j][i] = -1;
  }

  return 0;
}

int initialize_equilibrium(equilibrium_t *e)
{
  int i;

  /* verbose level one by default */
  e->verbose = 1;

  /* no equilibrium yet */
  e->isequil = 0;
  /* the state have not been set */
  e->is_state_set = 0;

  /* allocate a product structure */
  e->p = (product_t *)     malloc (sizeof(product_t));
  /* allocate a composition structure */
  e->c = (composition_t *) malloc (sizeof(composition_t));

  /* allocate the vector containing the delta ln(nj) */
  e->delta_ln_nj = (double *) calloc (200, sizeof(double));


  /* allocate the element vector and initialize to -1*/
  e->element = (int *) malloc (sizeof(int) * 100);
  for (i = 0; i < 100; i++)
    e->element[i] = -1;

  /* the composition have not been set */
  e->c->ncomp = 0;


  return 0;
}

int dealloc_equillibrium(equilibrium_t *e)
{
  
  dealloc_product (e->p);

  free (e->p);
  free (e->c);
  free (e->delta_ln_nj);
  free (e->element);

  free (e);

  return 0;
}

int dealloc_product(product_t *p)
{
  int i;
  for (i = 0; i < STATE_LAST; i++)
  {
    free (p->species[i]);
    free (p->coef[i]);
  }
  return 0;
}

int set_state(equilibrium_t *e, double T, double P)
{
  e->T = T;
  e->P = P;
  return 0;
}

int add_in_propellant(equilibrium_t *e, int sp, double mol)
{
  e->c->molecule[ e->c->ncomp ] = sp;
  e->c->coef[ e->c->ncomp ]     = mol;
  e->c->ncomp++;
  return 0;
}


int product_element_coef(int element, int molecule)
{
  int i;
  for (i = 0; i < 5; i++)
  {
    if (thermo_list[molecule].elem[i] == element)
      return thermo_list[molecule].coef[i];
  }
  return 0;
}

int propellant_element_coef(int element, int molecule)
{
  int i;
  for (i = 0; i < 6; i++)
  {
    if (propellant_list[molecule].elem[i] == element)
      return propellant_list[molecule].coef[i];
  }
  return 0;
}

int set_verbose(equilibrium_t *e, int v)
{
  e->verbose = v;
  return 0;
}


/* use the theory explain in
Theorical Evaluation of chemical propellant
by Roger Lawrence Wilkins
Prentice-Hall International series in space technology */
int fill_matrix(double **matrix, equilibrium_t *e)
{

  int i, j, k;
  double tmp, mol, tmp2;

  mol = 0.0;
  for (k = 0; k < e->p->n[GAS]; k++)
    mol += e->p->coef[GAS][k]; 


  // fill the matrix (part with the Lagrange multipliers)
  for (i = 0; i < e->n_element; i++) // each column
  {
    for (j = 0; j < e->n_element; j++) // each row
    {
      tmp = 0.0;
      for (k = 0; k < e->p->n[GAS]; k++) {
  	tmp += product_element_coef( e->element[j], e->p->species[GAS][k]) * 
	  product_element_coef( e->element[i], e->p->species[GAS][k] ) * 
	  e->p->coef[GAS][k];
      }
      matrix[j][i] = tmp;
    }
  }
  
  // X1c
  for (i = 0; i < e->p->n[CONDENSED]; i++) // column
  {
    for (j = 0; j < e->n_element; j++) // row
    {
      matrix[j][i + e->n_element ] = product_element_coef(e->element[j], 
					    e->p->species[CONDENSED][i]);
    }
  } 

  // delta ln(n)
  for (j = 0; j < e->n_element; j++)
  {
    tmp = 0.0;
    for (k = 0; k < e->p->n[GAS]; k++) {
  	tmp += product_element_coef( e->element[j], e->p->species[GAS][k]) * 
	  e->p->coef[GAS][k];
    }
    matrix[j][ e->n_element + e->p->n[CONDENSED] ] = tmp;
  }

  // right side
  for (j = 0; j < e->n_element; j++)
  {
    tmp = 0.0;

    for (k = 0; k < e->p->n[GAS]; k++)
      tmp += product_element_coef( e->element[j], e->p->species[GAS][k]) * 
	e->p->coef[GAS][k] * gibbs( e->p->species[GAS][k], GAS, 
				    e->p->coef[GAS][k],
				    e->n, e->T, e->P);

    // b[i]
    tmp2 =0.0;
    for (k = 0; k < STATE_LAST; k++)
      for (i = 0; i < e->p->n[k]; i++)
	tmp2 += product_element_coef( e->element[j], e->p->species[k][i]) *
	  e->p->coef[k][i];
 
    tmp -= tmp2;

    // b[i]o
    tmp2 = 0.0;
    for (i = 0; i < e->c->ncomp; i++)
      tmp2 += propellant_element_coef( e->element[j], e->c->molecule[i]) *
	e->c->coef[i];

    tmp += tmp2;

    matrix[j][ e->n_element + e->p->n[CONDENSED] + 1 ] = tmp;  
  }

  //second row
  for (i = 0; i < e->n_element; i++) // column
  {
    for (j = 0; j < e->p->n[CONDENSED]; j++) // row
    {
      //      printf("%d\n",j);
      matrix[j + e->n_element ][i] = product_element_coef(e->element[i],
					       e->p->species[CONDENSED][j]);
    }
  }
  
  /* set to zero */
  for (i = 0; i < e->p->n[CONDENSED]+1; i++) // column
  {
    for (j = 0; j < e->p->n[CONDENSED]; j++) // row
    {
      matrix[j + e->n_element ][i + e->n_element] = 0;
    }
  }

  // right side
  for (j = 0; j < e->p->n[CONDENSED]; j++) // row
  {
    matrix[ j + e->n_element ][ e->n_element + e->p->n[CONDENSED] + 1] =
      gibbs( e->p->species[CONDENSED][j], CONDENSED, e->p->coef[CONDENSED][j],
    	     e->n, e->T, e->P);  
  }

  // third big row
  for (i = 0; i < e->n_element; i++) // each column
  {   
      tmp = 0.0;
      for (k = 0; k < e->p->n[GAS]; k++)
  	tmp += product_element_coef( e->element[i], e->p->species[GAS][k] ) * 
	  e->p->coef[GAS][k];

      matrix[ e->n_element + e->p->n[CONDENSED] ][i] = tmp;      
  }

  // delta ln(n)
  matrix[e->n_element + e->p->n[CONDENSED]][e->n_element + e->p->n[CONDENSED]]
    =  mol - e->n;

  // right side
  tmp = 0.0;
  for (k = 0; k < e->p->n[GAS]; k++)
  {
    tmp += e->p->coef[GAS][k]*gibbs( e->p->species[GAS][k], GAS, 
				     e->p->coef[GAS][k],
				     e->n, e->T, e->P);
  }

  matrix[e->n_element + e->p->n[CONDENSED] ][e->n_element + e->p->n[CONDENSED] 
					  + 1] = e->n - mol + tmp;

  return 0;
}


/* may be optimize!!!!!!!! */
int remove_condensed(int *size, int *n, equilibrium_t *e)
{

  int i, j;
  int r = 0;

  for (i = 0; i < e->p->n[CONDENSED]; i++)
  {
    //printf("%s \n", (thermo_list + e->p->species[CONDENSED][i])->name );

    if (e->p->coef[CONDENSED][i] <= 0.0)
    {
      printf("%s should be remove\n", 
	     (thermo_list + e->p->species[CONDENSED][i])->name );
      /* remove from the list */
      for (j = i; j < e->p->n[CONDENSED]; j++)
      {
	e->p->species[CONDENSED][j] = 
	  e->p->species[CONDENSED][j + 1];
      }
      (*n)--;
      e->p->n[CONDENSED] = e->p->n[CONDENSED] - 1;
      (*size)--;
      r = 1;
    }
  }
  /* 0 if none remove */
  return r;
}

int include_condensed(equilibrium_t *e, int *n, int *size,
		      double **matrix, double *sol)
{
  double tmp;
  double temp;
  int i, j, k;
  int pos;

  tmp = 0.0;
  j   = -1;

  for (i = e->p->n[CONDENSED]; i < (*n); i++)
  {
    temp = 0.0;
    for (k = 0; k < e->n_element; k++)
      temp += sol[k]*product_element_coef(e->element[k], 
					  e->p->species[CONDENSED][i]);

    if ( gibbs0( e->p->species[CONDENSED][i], e->T) - temp < tmp )
    {
      tmp = gibbs0( e->p->species[CONDENSED][i], e->T) - temp;
      j = i; 
    }
  }

  if (!(j == -1))
  {
    printf("%s should be include\n", 
	   (thermo_list + e->p->species[CONDENSED][j])->name );

    /* to include the species, exchange the value */
    pos = e->p->species[CONDENSED][ e->p->n[CONDENSED] ];
    e->p->species[CONDENSED][ e->p->n[CONDENSED] ] = 
      e->p->species[CONDENSED][j];
    e->p->species[CONDENSED][j] = pos;
	
    e->p->n[CONDENSED]++;
    
    /* realloc the memory for matrix and sol */
    /* the size of the coefficient matrix */
    (*size) = e->n_element + e->p->n[CONDENSED] + 1;
  
    //printf("size: %d\n", *size);

    return 1;
  }

  return 0;
}

int equilibrium(equilibrium_t *equil)
{

  int       size;
  double ** matrix;
  double  * sol;
  int       convergence_ok;

  /* control factor */
  double lambda1;
  double lambda2; 
  double lambda;

  double mol; /* sommation of the gazeous nj */
  double temp;
  double tmplog;

  int i, j, k;

  int    n_condensed;


  print_propellant_composition(equil);

  list_element(equil);
  list_product(equil);

  n_condensed = equil->p->n[CONDENSED];

  /* initially, we do not consider the condensed */
  equil->p->n[CONDENSED] = 0;
  equil->n = 0.1; /* initial estimate of the mol number */

  /* the size of the coefficient matrix */
  size = equil->n_element + equil->p->n[CONDENSED] + 1;
  
  // allocate the memory for the matrix
  matrix = (double **) calloc (size, sizeof(double *));
  for (i = 0; i < size; i++)
    matrix[i] = (double *) calloc (size+1, sizeof(double));

  // allocate the memory for the solution vector
  sol = (double *) calloc (size, sizeof(double));



  mol = 0.0;
  for (i = 0; i < equil->p->n[GAS]; i++)
    mol += equil->p->coef[GAS][i];

  for (k = 0; k < iteration_max; k++)
  {

    fill_matrix(matrix, equil);

    if (equil->verbose > 1)
      print_matrix(matrix, size);

    if ( lu(matrix, sol, size) == -1)
    {
      printf("The matrix is singular, removing excess condensed.\n");
      /* remove undesired condensed */
      if (!remove_condensed(&size, &n_condensed, equil))
      {
	printf("None remove. Try reinserting remove gazes\n");
	for ( i = 0; i < equil->p->n[GAS]; i++)
	{
	  if (equil->p->coef[GAS][i] == 0.0)
	    equil->p->coef[GAS][i] = 1e-6;
	}
	k = -1;
	continue;	
      }
      //printf("size: %d\n", size);
      k = -1;
      continue;
    }


    if (equil->verbose > 1)
      print_vec(sol, size);

  
    /* compute the values of delta ln(nj) */
    equil->delta_ln_n = sol[ equil->n_element + equil->p->n[CONDENSED] ];

    for (i = 0; i < equil->p->n[GAS]; i++)
    {
      temp = 0.0;
      for (j = 0; j < equil->n_element; j++)
	temp += product_element_coef( equil->element[j], 
				      equil->p->species[GAS][i])*
	  sol[j];

   
      equil->delta_ln_nj[i] = -gibbs( equil->p->species[GAS][i], 
				      GAS, equil->p->coef[GAS][i], 
				      equil->n, equil->T, equil->P) 
	+ temp + equil->delta_ln_n;

      /* if the substance is exclude of the mixture */
      if ((equil->delta_ln_nj[i] < 0) && (equil->p->coef[GAS][i] == 0))
	equil->delta_ln_nj[i] = 0;

    }

    /* compute the control factor lambda */
    lambda1 = 0.0;
    for (i = 0; i < equil->p->n[GAS]; i++)
    {
      if ((equil->p->coef[GAS][i]/equil->n > conc_tol) && 
	  (equil->delta_ln_nj[i] > 0))
	lambda1 = _max( lambda1, fabs(equil->delta_ln_n), 
		       fabs(equil->delta_ln_nj[i]));
    }
    lambda1 = 2/lambda1;
    
    lambda2 = 1e12;
    for (i = 0; i < equil->p->n[GAS]; i++)
    {
      if ((equil->p->coef[GAS][i]/equil->n <= conc_tol) && 
	  (equil->delta_ln_nj[i]>0))
	lambda2 = __min(lambda2, fabs( ((-log(equil->p->coef[GAS][i]/equil->n) 
				       - 9.2103404)/(equil->delta_ln_nj[i] - 
						     equil->delta_ln_n))) );
    }
   
    lambda = _min(1, lambda1, lambda2);
    
  
    /* compute the new value for nj (gazeous) */
    
    if (equil->verbose > 1)
      printf(" \t  nj/n \t\t  Delta ln(nj)\n");
    for (i = 0; i < equil->p->n[GAS]; i++)
    {
      tmplog = log(equil->p->coef[GAS][i]) + lambda*equil->delta_ln_nj[i];
      
      if ( (tmplog - log(equil->n)) <= log(conc_tol))
	equil->p->coef[GAS][i] = 0.0;
      else
	equil->p->coef[GAS][i] = exp(tmplog);
      
      if (equil->verbose > 1)
	if (!(equil->p->coef[GAS][i] == 0))
	  printf("%s \t % .4e \t % .4e\n", 
		 (thermo_list + equil->p->species[GAS][i])->name, 
		 equil->p->coef[GAS][i]/equil->n,
		 equil->delta_ln_nj[i]);
    }
    

    /* compute the new value for nj (condensed) */
    for (i = 0; i < equil->p->n[CONDENSED]; i++)
    {
      equil->p->coef[CONDENSED][i] = equil->p->coef[CONDENSED][i] + 
	lambda*sol[equil->n_element + i];
      /* molar of condensed species */
      if (equil->verbose > 1)
	printf("%s: \t %f\n", 
	       (thermo_list + equil->p->species[CONDENSED][i])->name, 
	       equil->p->coef[CONDENSED][i]);
    }
    
    /* new value of n */
    equil->n = exp( log(equil->n) + lambda*equil->delta_ln_n );
    
    
    mol = 0.0;
    for (i = 0; i < equil->p->n[GAS]; i++)
      mol += equil->p->coef[GAS][i];

    
    // check for convergence 
    for (i = 0; i < equil->p->n[GAS]; i++)
    {
      if (!(equil->p->coef[GAS][i]*fabs(equil->delta_ln_nj[i])/mol <=conv_tol))
	convergence_ok = 0;
      
    }
    
    for ( i = 0; i < equil->p->n[CONDENSED]; i++ )
    {
    /* test for the condensed phase */
      if (!(sol[equil->n_element+1]/mol <= conv_tol))
	convergence_ok = 0;
    }

    if (!(fabs(equil->delta_ln_n) <= conv_tol))
      convergence_ok = 0;
    

    if (convergence_ok)
    {
      printf("The solution has converge in %d iteration\n\n", k);
      //print_product_composition(equil);

      /* remove undesired condensed */
      if (remove_condensed(&size, &n_condensed, equil))
      	convergence_ok = 0;
      
      /* find if a new condensed species should be include */
      if ( include_condensed(equil, &n_condensed, &size, matrix, sol) )
      {
	free(matrix);
	free(sol);
	// allocate the memory for the matrix
	matrix = (double **) malloc (sizeof(double *) * size);
	for (i = 0; i < size; i++)
	  matrix[i] = (double *) malloc (sizeof(double) * (size+1));
	
	// allocate the memory for the solution vector
	sol = (double *) malloc (sizeof(double) * size);
	
	k = -1;
	continue; /* begin a new loop */
      }

      if (convergence_ok)
	break;
      k = -1;
      continue;
    }
    else
      if (equil->verbose > 1)
	{
	  printf("The solution doesn't converge\n\n");
	  //remove_condensed(&size, &n_condensed, equil);
	}


    /* suppose that it will converge the next time */
    convergence_ok = 1;
    
  }

  print_product_composition(equil);

  free(sol);
  free (matrix);
  return 0;
}
  




