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
#include "type.h"

/* Initial temperature estimate for problem with not-fixed temperature */
#define ESTIMATED_T 3000

/* global variable containing the information about chemical species */
propellant_t	*propellant_list;
thermo_t	    *thermo_list;


double conc_tol   = 1e-15;
double conv_tol   = 0.5e-5;
int iteration_max = 200;

/* 1 for verbose, 0 for non-verbose */
int global_verbose = 0;

int global_error = 0;

/****************************************************************
VARIABLE: Contain the molar mass of element by atomic number
          molar_mass[0] contain hydrogen and so on.
          Data come from Sargent-Welch 1996
*****************************************************************/
const float molar_mass[N_SYMB] = { 
  1.00794,  4.00260, 6.941,    9.01218,  10.811,   12.011,
  14.0067,  15.9994, 18.99840, 20.11797, 22.99977, 24.305, 
  26.98154, 28.0855, 30.97376, 32.066,   35.4527,  39.948, 
  39.0983,  40.078,  44.9559,  47.88,    50.9415,  51.996,
  54.938,   55.847,  58.9332,  58.6934,  63.546,   65.39,
  69.723,   72.61,   74.9216,  78.96,    79.904,   83.80,
  85.4678,  87.62,   88.9059,  91.224,   92.9064,  95.94,
  98.0,     101.07,  102.9055, 106.42,   107.868,  112.41,
  114.82,   118.71,  121.757,  127.60,   126.9045, 131.29,
  132.9054, 137.33,  138.9055, 140.12,   140.9077, 144.24,
  145.,     150.36,  151.965,  157.25,   158.9253, 162.50,
  164.9303, 167.26,  168.9342, 173.04,   174.967,  178.49,
  180.9479, 183.85,  186.207,  190.2,    192.22,   195.08,
  196.9665, 200.59,  204.383,  207.2,    208.9804, 209.,
  210.,     222.,    223.,     226.0254, 227.,     232.0381,
  231.0359, 238.029, 237.0482, 244.,     12.011,   9.01218,
  10.811,   24.305,  26.98154, 257.0,    0,        2};


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
VARIABLE: molar gaz constant in J/(mol K)
****************************************************************/
const float  R = 8.3143; 

/* Enthalpy in the standard state */
double enthalpy_0(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);
  double val;
  int i, pos;

  if (T < s->range[0][0]) /* Temperature below the lower range */
  {
    pos = 0;
  }       /*Temperature above the higher range */
  else if (T >= s->range[s->nint-1][1]) 
  {
    pos = s->nint - 1;
  }
  else
  {
    for (i = 0; i < s->nint; i++) /* Find the range */
    {
      if ((T >= s->range[i][0]) && (T < s->range[i][1]))
        pos = i;
    }
  }
  
  /* parametric equation for dimentionless enthalpy */
  val = -s->param[pos][0]*pow(T, -2) + s->param[pos][1]*pow(T, -1)*log(T) +
    s->param[pos][2] + s->param[pos][3]*T/2 + s->param[pos][4]*pow(T, 2)/3 +
    s->param[pos][5]*pow(T, 3)/4 + s->param[pos][6]*pow(T, 4)/5 +
    s->param[pos][7]/T;

  return val; /* dimensionless enthalpy */
  
  //if (global_verbose > 1)
  //  printf("Error: temperature(%.0f) out of range for %d: %s\n", T,
  //         sp, s.name);
  //global_error = 1;
  //return -1;
}

/* Entropy in the stabdard state */
double entropy_0(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);
  double val;
  int i, pos;

  if (T < s->range[0][0])
  {
    pos = 0;
  }
  else if (T >= s->range[s->nint-1][1])
  {
    pos = s->nint - 1;
  }
  else
  {
    for (i = 0; i < s->nint; i++)
    {
      if ((T >= s->range[i][0]) && (T < s->range[i][1]))
        pos = i;
    }
  }
  
  /* parametric equation for dimentionless entropy */
  val = -s->param[pos][0]*pow(T, -2)/2 - s->param[pos][1]*pow(T, -1)
    + s->param[pos][2]*log(T) + s->param[pos][3]*T
    + s->param[pos][4]*pow(T, 2)/2
    + s->param[pos][5]*pow(T, 3)/6 + s->param[pos][6]*pow(T, 4)/4 
    + s->param[pos][8];
  
  return val;
  
  //if (global_verbose > 1)
  //  printf("Error: temperature(%.0f) out of range for %d: %s\n", T,
  //         sp, s.name);
  //global_error = 1;
  //return -1;
}

/* Specific heat in the standard state */
double specific_heat_0(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);
  double val;
  int i, pos;

  if (T < s->range[0][0])
  {
    pos = 0;
  }
  else if (T >= s->range[s->nint-1][1])
  {
    pos = s->nint - 1;
  }
  else
  {
    for (i = 0; i < s->nint; i++)
    {
      if ((T >= s->range[i][0]) && (T < s->range[i][1]))
        pos = i;
    }
  }
  
  /* parametric equation for dimentionless specific_heat */
  val = s->param[pos][0]*pow(T, -2) + s->param[pos][1]*pow(T, -1)
    + s->param[pos][2] + s->param[pos][3]*T + s->param[pos][4]*pow(T, 2)
    + s->param[pos][5]*pow(T, 3) + s->param[pos][6]*pow(T, 4);

  return val;
  
  // if (global_verbose > 1)
  // printf("Error: temperature(%.0f) out of range for %d: %s\n", T,
  //        sp, s.name);
  //global_error = 1;
  //return -1;
}

/* 0 if out of range, 1 if ok */
int temperature_check(int sp, float T)
{
  int i;
  thermo_t *s = (thermo_list + sp);
  for (i = 0; i < s->nint; i++)
  {
    if ((T >= s->range[i][0]) && (T < s->range[i][1]))
      return 1;
  }
  return 0;
}

/* unuse function */
double delta_enthalpy(int sp, float T)
{
  /* delta henthalpy in J/mol */
  return enthalpy_0(sp, T)*R*T - (thermo_list + sp)->heat;
}


double entropy(int sp, state_t st, double nj, double n, float T, float P)
{
  double s;
  
  if ( (st == GAS) && (nj == 0.0))
    return 0;
  
  switch (st)
  {
    case GAS:    
        s = entropy_0(sp, T) - log(nj/n) - log(P);
        break;
    case CONDENSED:
        s = entropy_0(sp, T);
        break;
    default:
        s = 0;
  }  
  return s;
}


/* uo = HO(T) - S(T)  (from the nasa book) */
double gibbs_0(int sp, float T)
{
  return enthalpy_0(sp, T) - entropy_0(sp, T); /* dimensionless */
}

/* J/mol T is in K, P is in atm */
double gibbs(int sp, state_t st, double nj, double n, float T, float P)
{
  double g;

  if ( (st == GAS) && (nj == 0.0))
    return 0;

  switch (st)
  {
    case GAS:    
        g = gibbs_0(sp, T) + log(nj/n) + log(P);
        break;
    case CONDENSED:
        g = gibbs_0(sp, T);
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
    ans += (propellant_list + molecule)->coef[i] *
      molar_mass[(propellant_list + molecule)->elem[i]];
  return ans;
}

/* J/mol */
double heat_of_formation(int molecule)
{
  double hf = (propellant_list + molecule)->heat * 
    propellant_molar_mass(molecule);
  return hf;
}


double propellant_enthalpy(equilibrium_t *e) //composition_t *c)
{
  int i;
  double h = 0.0;
  for (i = 0; i < e->c->ncomp; i++)
  {
    h += e->c->coef[i]*heat_of_formation( e->c->molecule[i] )/(R*e->T);
  }
  return h;
}

double product_enthalpy(equilibrium_t *e)
{
  int i;
  double h = 0.0;
  /* for gases */
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    h += e->p->coef[GAS][i]*enthalpy_0(e->p->species[GAS][i], e->T);
  }
  /* for condensed */
  for (i = 0; i < e->p->n[CONDENSED]; i++)
  {
    h += e->p->coef[CONDENSED][i]*enthalpy_0(e->p->species[CONDENSED][i], e->T);
  }
  return h;
}

double product_entropy(equilibrium_t *e)
{
  int i;
  double ent = 0.0;
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    ent += e->p->coef[GAS][i]*entropy(e->p->species[GAS][i], GAS,
                                      e->p->coef[GAS][i],
                                      e->n, e->T, e->P);
  }
  for (i = 0; i < e->p->n[CONDENSED]; i++)
  {
    ent += e->p->coef[CONDENSED][i]*entropy(e->p->species[CONDENSED][i],
                                            CONDENSED,
                                            e->p->coef[CONDENSED][i],
                                            e->n, e->T, e->P);
  }
  return ent;
}

int thermo_search(char *str)
{
  int i;
  int last = -1;
  
  for (i = 0; i < MAX_THERMO; i++)
  {
    if (!(strncasecmp(str, (thermo_list + i)->name, strlen(str))))
      {
        last = i;
        printf("%-5d %s\n", i, (thermo_list + i)->name);
      }
  }
  return last;
}

int propellant_search(char *str)
{
  int i;
  int last = -1;
  
  for (i = 0; i < MAX_PROPELLANT; i++)
  {
    if (!(strncasecmp(str, (propellant_list + i)->name, strlen(str))))
    {
      last = i;
      printf("%-5d %s\n", i, (propellant_list + i)->name);
    }
  }
  return last;
  
}

double product_molar_mass(equilibrium_t *e)
{
  int i;
  double mol  = 0.0;
  double mass = 0.0;

  for (i = 0; i < e->p->n[GAS]; i++)
  {
    mol  += e->p->coef[GAS][i]; 
    mass += e->p->coef[GAS][i]*(thermo_list + e->p->species[GAS][i])->weight;  
  }
  
  for (i = 0; i < e->p->n[CONDENSED]; i++)
  {
    mass += e->p->coef[CONDENSED][i]*(thermo_list + 
                                      e->p->species[CONDENSED][i])->weight;
  }
  return (mass/mol);
}

int list_element(equilibrium_t *e)
{
  int n = 0;
  int t = 0;
  int i, j, k;
  
  for (i = 0; i < e->c->ncomp; i++)
  {
    /* maximum of 6 different atoms in the composition */
    for (j = 0; j < 6; j++)
    {	       
      if (!( (propellant_list + e->c->molecule[i])->coef[j] == 0))
      {
        /* get the element */
        t = (propellant_list + e->c->molecule[i])->elem[j];
        
        for (k = 0; k <= n; k++)
        {
          /* verify if the element was not already in the list */
          if (e->element[k] == t)
            break;
          /* if we have check each element, add it to the list */
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

  if (e->verbose > 0)
  {
    printf("%d different elements in the propellant\n", n);
    /* Print those elements */
    for (i = 0; i < n; i++)
      printf("%s ", symb[e->element[i]] );
    printf("\n");
  }
  
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
    /* for each of the five possible element of a species */
    for (k = 0; k < 5; k++)
    {
      if (!((thermo_list + j)->coef[k] == 0))
      {
        for (i = 0; i < e->n_element; i++)
        {
          if (e->element[i] == (thermo_list + j)->elem[k])
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
      st = (thermo_list + j)->state;

      /* maybe useful to add every possible product  || 1 */ 
      if (temperature_check(j, e->T))  // if the molecule could 
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
    if ((e->p->species[i] = (int *)realloc (e->p->species[i], 
                                            sizeof(int) * e->p->n[i]))
        == NULL)
      printf("Reallocation of memory failed\n");
    
    if ((e->p->coef[i] = (double *)realloc (e->p->coef[i], 
                                            sizeof(double) * e->p->n[i]))
        == NULL)
      printf("Reallocation of memory failed\n");
  }
  
  /* initialize tho mol number to 0.1mol/(nb of gazeous species) */
  for (i = 0; i < e->p->n[GAS]; i++)
    e->p->coef[GAS][i] = 0.1/e->p->n[GAS];
  
  /* initialize condensed to zero */
  for (i = 0; i < e->p->n[CONDENSED]; i++)
    e->p->coef[CONDENSED][i] = 0;

  if (e->verbose > 0)
  {
    printf("%d possible combustion product\n", n);
    printf("%d gazeous species\n", e->p->n[GAS]);
    if (e->verbose > 1)
      print_gazeous(*(e->p));
    printf("%d condensed species\n", e->p->n[CONDENSED]);
    if (e->verbose > 1)
      print_condensed(*(e->p));
  }
  
  return n;
}


int mem_alloc(void)
{
  if ((thermo_list = 
       (thermo_t *)malloc (sizeof(thermo_t) * MAX_THERMO)) == NULL)
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
  // the maximum of 300 species could be not enough
  for (i = 0; i < STATE_LAST; i++)
    p->species[i] = (int *)malloc (sizeof(int) * 300);
  
  for (i = 0; i < STATE_LAST; i++)
    p->coef[i] = (double *)malloc (sizeof(double) * 300);
  
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
  e->delta_ln_nj = (double *)calloc (200, sizeof(double));

  /* allocate the element vector and initialize to -1*/
  e->element = (int *)malloc (sizeof(int) * 100);
  for (i = 0; i < 100; i++)
    e->element[i] = -1;

  /* the composition have not been set */
  e->c->ncomp = 0;

  /* initialize the product */
  initialize_product(e->p);

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

/* Mass of propellant in gram */
double propellant_mass(equilibrium_t *e)
{
  int i;
  double mass = 0.0;
  for (i = 0; i < e->c->ncomp; i++)
  {
    mass += e->c->coef[i]*propellant_molar_mass(e->c->molecule[i]);
  }
  return mass;
}

int product_element_coef(int element, int molecule)
{
  int i;
  for (i = 0; i < 5; i++)
  {
    if ((thermo_list + molecule)->elem[i] == element)
      return (thermo_list + molecule)->coef[i];
  }
  return 0;
}

int propellant_element_coef(int element, int molecule)
{
  int i;
  for (i = 0; i < 6; i++)
  {
    if ((propellant_list + molecule)->elem[i] == element)
      return (propellant_list + molecule)->coef[i];
  }
  return 0;
}

int set_verbose(equilibrium_t *e, int v)
{
  e->verbose = v;
  return 0;
}


/* Compute an initial estimate of the product composition using
   a method develop by G. Eriksson */
int initial_estimate(equilibrium_t *e)
{
  int i, j, mol;
  int components = e->p->n[GAS] + e->p->n[CONDENSED];
  double energy[components];

  for (i = 0; i < e->p->n[CONDENSED]; i++)
  {
    mol = 0;
    for (j = 0; j < 5; j++)
      mol += (thermo_list + e->p->species[CONDENSED][i])->coef[j];

    energy[i] = (enthalpy_0(e->p->species[CONDENSED][i], e->T) -
                 entropy_0(e->p->species[CONDENSED][i], e->T))*R*e->T/mol;
    printf("%s \t %f \t %i\n",
           (thermo_list + e->p->species[CONDENSED][i])->name,
           energy[i], mol);
  }
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    mol = 0;
    for (j = 0; j < 5; j++)
      mol += (thermo_list + e->p->species[GAS][i])->coef[j];
    
    energy[i + e->p->n[CONDENSED]] = (enthalpy_0(e->p->species[GAS][i],
                                               e->T) -
                                      entropy_0(e->p->species[GAS][i],
                                              e->T))*R*e->T/mol;
    
    printf("%s \t %f \t %i\n",
           (thermo_list + e->p->species[GAS][i])->name,
           energy[i + e->p->n[CONDENSED]], mol);
  }
  
  return 0;
}

/* use the theory explain in
Theorical Evaluation of chemical propellant
by Roger Lawrence Wilkins
Prentice-Hall International series in space technology */
int fill_matrix(double **matrix, equilibrium_t *e, problem_t P)
{

  int i, j, k;
  double tmp, mol;

  /* position of the right side dependeing on the type of problem */
  int roff = 2;
 
  if (P == TP)
    roff = 1;

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
  
  // Delta n
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
  
  // delta ln(T) (for SP and HP only)
  if (P != TP)
  {
    for (j = 0; j < e->n_element; j++)
    {
      tmp = 0.0;
      for (k = 0; k < e->p->n[GAS]; k++) {
        tmp += product_element_coef( e->element[j], e->p->species[GAS][k]) *
          e->p->coef[GAS][k] * enthalpy_0( e->p->species[GAS][k], e->T);
      }
      matrix[j][ e->n_element + e->p->n[CONDENSED] + 1] = tmp;
    }
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
    for (k = 0; k < STATE_LAST; k++)
      for (i = 0; i < e->p->n[k]; i++)
        tmp -= product_element_coef( e->element[j], e->p->species[k][i]) *
          e->p->coef[k][i];
    
    // b[i]o
    for (i = 0; i < e->c->ncomp; i++)
      tmp += propellant_element_coef( e->element[j], e->c->molecule[i]) *
        e->c->coef[i];
    
    matrix[j][ e->n_element + e->p->n[CONDENSED] + roff ] = tmp;  
  }
  
  //second row
  for (i = 0; i < e->n_element; i++) // column
  {
    for (j = 0; j < e->p->n[CONDENSED]; j++) // row
    {
      /* copy the symetric part of the matrix */
      matrix[j + e->n_element ][i] = matrix[i][j + e->n_element ];
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
  
  /* delta ln(T) */
  if (P != TP)
  {
    for (j = 0; j < e->p->n[CONDENSED]; j++) // row
      matrix[ j + e->n_element ][ e->n_element + e->p->n[CONDENSED] + 1] = 
        enthalpy_0( e->p->species[CONDENSED][j], e->T);
  }


  // right side
  for (j = 0; j < e->p->n[CONDENSED]; j++) // row
  {
    matrix[ j + e->n_element ][ e->n_element + e->p->n[CONDENSED] + roff] =
      gibbs( e->p->species[CONDENSED][j], CONDENSED, e->p->coef[CONDENSED][j],
             e->n, e->T, e->P);  
  }

  // third big row
  for (i = 0; i < e->n_element; i++) // each column
  {   
    /* copy the symetric part of the matrix */
    matrix[ e->n_element + e->p->n[CONDENSED] ][i] =
      matrix[i][ e->n_element + e->p->n[CONDENSED] ];
  }

  /* set to zero */
  for (i = 0; i < e->p->n[CONDENSED]; i++) // column
  {
      matrix[e->n_element + e->p->n[CONDENSED] ][i + e->n_element] = 0;
  }
  
  // delta ln(n)
  matrix[e->n_element + e->p->n[CONDENSED]][e->n_element + e->p->n[CONDENSED]]
    =  mol - e->n;

  /* delta ln(T) */
  tmp = 0.0;
  for (k = 0; k < e->p->n[GAS]; k++)
    tmp += e->p->coef[GAS][k]*enthalpy_0( e->p->species[GAS][k], e->T ); 
  matrix[e->n_element + e->p->n[CONDENSED]][e->n_element + e->p->n[CONDENSED] 
                                           + 1] = tmp;

  // right side
  tmp = 0.0;
  for (k = 0; k < e->p->n[GAS]; k++)
  {
    tmp += e->p->coef[GAS][k]*gibbs( e->p->species[GAS][k], GAS, 
                                     e->p->coef[GAS][k],
                                     e->n, e->T, e->P);
  }

  matrix[e->n_element + e->p->n[CONDENSED] ][e->n_element + e->p->n[CONDENSED] 
                                            + roff] = e->n - mol + tmp;

  // for enthalpy/pressure problem
  if (P == HP)
  {
    /* part with lagrangian multipliers */
    for (i = 0; i < e->n_element; i++) // each column
    {   
      tmp = 0.0;
      for (k = 0; k < e->p->n[GAS]; k++)
        tmp += product_element_coef( e->element[i], e->p->species[GAS][k] ) * 
          e->p->coef[GAS][k] * enthalpy_0( e->p->species[GAS][k], e->T);
      
      matrix[ e->n_element + e->p->n[CONDENSED] + 1 ][i] = tmp;      
    }

    // Delta n
    for (i = 0; i < e->p->n[CONDENSED]; i++)
      matrix[ e->n_element + e->p->n[CONDENSED] + 1 ][i + e->n_element] = 
        enthalpy_0( e->p->species[CONDENSED][i], e->T);

    // Delta ln(n)
    tmp = 0.0;
    for (k = 0; k < e->p->n[GAS]; k++)
      tmp += e->p->coef[GAS][k]*enthalpy_0( e->p->species[GAS][k], e->T ); 

    matrix[e->n_element + e->p->n[CONDENSED] + 1][e->n_element + 
                                                 e->p->n[CONDENSED] ] = tmp;
    

    tmp = 0.0;
    for (k = 0; k < e->p->n[GAS]; k++)
      tmp += e->p->coef[GAS][k]*specific_heat_0( e->p->species[GAS][k], e->T );

    for (k = 0; k < e->p->n[CONDENSED]; k++)
      tmp += e->p->coef[CONDENSED][k]*
        specific_heat_0( e->p->species[CONDENSED][k], e->T);

    for (k = 0; k < e->p->n[GAS]; k++)
      tmp += e->p->coef[GAS][k]*enthalpy_0( e->p->species[GAS][k], e->T)*
        enthalpy_0( e->p->species[GAS][k], e->T);
    
    matrix[e->n_element + e->p->n[CONDENSED] + 1][e->n_element +
                                                 e->p->n[CONDENSED] + 1] = tmp;

    /* right side */
    tmp = 0.0;
    /* enthalpy of reactant (not sure for the value of R */
    //for (k = 0; k < e->c->ncomp; k++)
    //  tmp += e->c->coef[k]*heat_of_formation( e->c->molecule[k] );
 
    //for (k = 0; k < e->p->n[GAS]; k++)
    //  tmp -= e->p->coef[GAS][k]*enthalpy( e->p->species[GAS][k], e->T);

    //for (k = 0; k < e->p->n[CONDENSED]; k++)
    //  tmp -= e->p->coef[CONDENSED][k]*enthalpy( e->p->species[CONDENSED][k], 
    //                                            e->T);

    tmp = propellant_enthalpy(e) - product_enthalpy(e);
    
    for (k = 0; k < e->p->n[GAS]; k++)
      tmp += e->p->coef[GAS][k]*enthalpy_0( e->p->species[GAS][k], e->T)*
        gibbs( e->p->species[GAS][k], GAS, e->p->coef[GAS][k], e->n, 
               e->T, e->P);
    
    matrix[e->n_element + e->p->n[CONDENSED] + 1][e->n_element +
                                                 e->p->n[CONDENSED] + 2] = tmp;
    
  } // for entropy/pressure problem
  else if (P == SP)
  {
    /* part with lagrangian multipliers */
    for (i = 0; i < e->n_element; i++) // each column
    {   
      tmp = 0.0;
      for (k = 0; k < e->p->n[GAS]; k++)
        tmp += product_element_coef( e->element[i], e->p->species[GAS][k] ) * 
          e->p->coef[GAS][k] * entropy(e->p->species[GAS][i], GAS,
                                       e->p->coef[GAS][i],
                                       e->n, e->T, e->P);
      
      matrix[ e->n_element + e->p->n[CONDENSED] + 1][i] = tmp;      
    }
    
    // Delta n
    for (i = 0; i < e->p->n[CONDENSED]; i++)
      matrix[ e->n_element + e->p->n[CONDENSED] + 1 ][i + e->n_element] = 
        entropy_0( e->p->species[CONDENSED][i], e->T); /* ok for condensed */
    

    // Delta ln(n)
    tmp = 0.0;
    for (k = 0; k < e->p->n[GAS]; k++)
      tmp += e->p->coef[GAS][k]*entropy(e->p->species[GAS][i], GAS,
                                        e->p->coef[GAS][i],
                                        e->n, e->T, e->P);

    matrix[e->n_element + e->p->n[CONDENSED] + 1][e->n_element + 
                                                 e->p->n[CONDENSED] ] = tmp;
    
    tmp = 0.0;
    for (k = 0; k < e->p->n[GAS]; k++)
      tmp += e->p->coef[GAS][k]*specific_heat_0( e->p->species[GAS][k], e->T );

    for (k = 0; k < e->p->n[CONDENSED]; k++)
      tmp += e->p->coef[CONDENSED][k]*
        specific_heat_0( e->p->species[CONDENSED][k], e->T);

    for (k = 0; k < e->p->n[GAS]; k++)
      tmp += e->p->coef[GAS][k]*enthalpy_0( e->p->species[GAS][k], e->T)*
        entropy(e->p->species[GAS][i], GAS,
                e->p->coef[GAS][i],
                e->n, e->T, e->P);
    
    matrix[e->n_element + e->p->n[CONDENSED] + 1][e->n_element +
                                                 e->p->n[CONDENSED] + 1] = tmp;
    
    
    tmp = 0.0;
    /* entropy of reactant (not sure for the value of R */
    tmp += 132; // assign entropy

    tmp -= product_entropy(e);
      
      //for (k = 0; k < e->p->n[GAS]; k++)
      //tmp -= e->p->coef[GAS][k]*entropy_0( e->p->species[GAS][k], e->T);

    tmp += e->n;

    for (k = 0; k < e->p->n[GAS]; k++)
      tmp -= e->p->coef[GAS][k];

    for (k = 0; k < e->p->n[GAS]; k++)
      tmp += e->p->coef[GAS][k]*gibbs(e->p->species[GAS][k], GAS,
                                      e->p->coef[GAS][k], e->n, e->T, e->P)*
        entropy(e->p->species[GAS][i], GAS,
                e->p->coef[GAS][i],
                e->n, e->T, e->P);
    
    matrix[e->n_element + e->p->n[CONDENSED] + 1][e->n_element +
                                                 e->p->n[CONDENSED] + 2] = tmp;
    
  }
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

    if (e->p->coef[CONDENSED][i] <= 0.0 ||
        !(temperature_check(e->p->species[CONDENSED][i], e->T)) )
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
    if (temperature_check(e->p->species[CONDENSED][i], e->T))
    {
      temp = 0.0;
      for (k = 0; k < e->n_element; k++)
        temp += sol[k]*product_element_coef(e->element[k], 
                                            e->p->species[CONDENSED][i]);
      
      if ( gibbs_0( e->p->species[CONDENSED][i], e->T) - temp < tmp )
      {
        tmp = gibbs_0( e->p->species[CONDENSED][i], e->T) - temp;
        j = i; 
      }
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
    
    return 1;
  }
  return 0;
}


int new_approximation(equilibrium_t *e, double *sol, problem_t P)
{
  
  int i, j;

  /* control factor */
  double lambda1, lambda2, lambda;

  double temp;
  
  /* compute the values of delta ln(nj) */
  e->delta_ln_n = sol[ e->n_element + e->p->n[CONDENSED] ];
      
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    temp = 0.0;
    for (j = 0; j < e->n_element; j++)
      temp += product_element_coef(e->element[j], 
                                   e->p->species[GAS][i]) * sol[j];
        
    e->delta_ln_nj[i] = -gibbs(e->p->species[GAS][i], 
                               GAS, e->p->coef[GAS][i], 
                               e->n, e->T, e->P) 
      + temp + e->delta_ln_n;
        
    /* if the substance is exclude of the mixture */
    if ((e->delta_ln_nj[i] < 0) && (e->p->coef[GAS][i] == 0))
      e->delta_ln_nj[i] = 0;
        
  }
      
  /* compute the control factor lambda */
  lambda1 = 0.0;
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    if ((e->p->coef[GAS][i]/e->n > conc_tol) && 
        (e->delta_ln_nj[i] > 0))
      lambda1 = _max( lambda1, fabs(e->delta_ln_n), 
                      fabs(e->delta_ln_nj[i]));
  }
  lambda1 = 2/lambda1;
      
  lambda2 = 1e12; /* a big number */
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    if ((e->p->coef[GAS][i]/e->n <= conc_tol) && 
        (e->delta_ln_nj[i]>0))
      lambda2 = __min(lambda2, fabs( ((-log(e->p->coef[GAS][i]/e->n) 
                                       - 9.2103404)/(e->delta_ln_nj[i] - 
                                                     e->delta_ln_n))) );
  }
  
  lambda = _min(1, lambda1, lambda2);
         
  /* compute the new value for nj (gazeous) */
      
  if (e->verbose > 1)
    printf(" \t  nj/n \t\t  Delta ln(nj)\n");
      
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    temp = log(e->p->coef[GAS][i]) + lambda*e->delta_ln_nj[i];
    
    if ( (temp - log(e->n)) <= log(conc_tol))
      e->p->coef[GAS][i] = 0.0;
    else
      e->p->coef[GAS][i] = exp(temp);
        
    if (e->verbose > 1)
      if (!(e->p->coef[GAS][i] == 0))
        printf("%s \t % .4e \t % .4e\n", 
               (thermo_list + e->p->species[GAS][i])->name, 
               e->p->coef[GAS][i]/e->n,
               e->delta_ln_nj[i]);
  }
          
  /* compute the new value for nj (condensed) */
  for (i = 0; i < e->p->n[CONDENSED]; i++)
  {
    e->p->coef[CONDENSED][i] = e->p->coef[CONDENSED][i] + 
      lambda*sol[e->n_element + i];
        
    /* molar of condensed species */
    if (e->verbose > 1)
      printf("%s: \t %f\n", 
             (thermo_list + e->p->species[CONDENSED][i])->name, 
             e->p->coef[CONDENSED][i]);
  }
  
  /* new value of T */
  if (P != TP)
    e->T = exp( log(e->T) + lambda*sol[e->n_element + 
                                      e->p->n[CONDENSED] + 1]);
      
  if (e->verbose > 1)
    printf("Temperature: %f\n", e->T);
      
  /* new value of n */
  e->n = exp( log(e->n) + lambda*e->delta_ln_n );
      
  return 0;
}
    
bool convergence(equilibrium_t *e, double *sol)
{
  int i;
  double mol;
  
  /* compute the total number of mol */
  mol = 0.0;
  for (i = 0; i < e->p->n[GAS]; i++)
    mol += e->p->coef[GAS][i];
      
      
  // check for convergence 
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    if (!(e->p->coef[GAS][i]*fabs(e->delta_ln_nj[i])/mol <= conv_tol))
      return false; /* haven't converge yet */
  }
      
  for ( i = 0; i < e->p->n[CONDENSED]; i++ )
  {
    /* test for the condensed phase */
    if (!(sol[e->n_element+1]/mol <= conv_tol))
      return false; /* haven't converge yet */
  }
      
  if (!(fabs(e->delta_ln_n) <= conv_tol))
    return false; /* haven't converge yet */

  return true;
}

int equilibrium(equilibrium_t *equil, problem_t P)
{

  int       size;     /* size of the matrix */
  double ** matrix;   
  double  * sol;

  bool      convergence_ok;
  bool      stop           = false;
  bool      gas_reinserted = false;
  bool      solution_ok    = false;

  int i, k;

  int n_condensed;  /* number of condensed species */

  /* position of the right side dependeing on the type of problem */
  int roff = 2;
  
  if (P == TP)
    roff = 1;

  /* initial temperature for assign enthalpy, entropy/pressure */
  if ( P != TP)
    equil->T = ESTIMATED_T;
 
  print_propellant_composition(equil);

  list_element(equil);
  list_product(equil);

  printf("----------------------\n");
  printf("Beginning computation\n");

  n_condensed = equil->p->n[CONDENSED];

  /* First determine an initial estimate of the composition
     to accelerate the convergence */
  /* initial_estimate(equil); */
  
  /* initially, we do not consider the condensed */
  equil->p->n[CONDENSED] = 0;
  equil->n = 0.1; /* initial estimate of the mol number */

  /* the size of the coefficient matrix */
  size = equil->n_element + equil->p->n[CONDENSED] + roff;
  
  // allocate the memory for the matrix
  matrix = (double **) malloc (sizeof(double *) * size);
  for (i = 0; i < size; i++)
    matrix[i] = (double *) malloc (sizeof(double) * (size+1));
  
  // allocate the memory for the solution vector
  sol = (double *) calloc (size, sizeof(double));

 
  /* main loop */
  for (k = 0; k < iteration_max; k++)
  {

    /* Initially we haven't a good solution */
    solution_ok = false;

    while (!solution_ok)
    {
      global_error = 0;
      fill_matrix(matrix, equil, P);

      if (global_error)
        printf("Global error occur\n"); /* probably a condensed that could
                                           not exist at the temperature */
      
      if (equil->verbose > 1)
        print_matrix(matrix, size);
    
      if ( lu(matrix, sol, size) == -1)
      {
        /* the matrix have no unique solution */
        printf("The matrix is singular, removing excess condensed.\n");
          
        /* Try removing excess condensed */
        if (!remove_condensed(&size, &n_condensed, equil))
        {
          if (gas_reinserted)
          {
            printf("No convergence, don't trust results\n");
            /* finish the main loop */
            stop = true;
            break;
            
          }
          printf("None remove. Try reinserting remove gaz\n");
          for ( i = 0; i < equil->p->n[GAS]; i++)
          {
            if (equil->p->coef[GAS][i] == 0.0)
              equil->p->coef[GAS][i] = 1e-6;
          }
          gas_reinserted = true;
        }
        else
          gas_reinserted = false;
          
          /* Restart the loop counter to zero for a new loop */
        k = 0;
      }
      else
      {
        /* There is a solution to the matrix */
        solution_ok = true;
      }
        
    }
      
    if (equil->verbose > 1)
      print_vec(sol, size);    /* print the solution vector */

    /* compute the new approximation */
    new_approximation(equil, sol, P);
    
    convergence_ok = false;

    /* verify the convergence */
    
    if (convergence(equil, sol))
    {
      convergence_ok = true;
      
      printf("The solution has converge in %d iteration\n", k);
        
      gas_reinserted = false;
        
      /* remove undesired condensed */
      if (remove_condensed(&size, &n_condensed, equil))
        convergence_ok = false;  /* condensed have been removed */
        
      /* find if a new condensed species should be include */
      if (include_condensed(equil, &n_condensed, &size, matrix, sol))
      {
        free(matrix);
        free(sol);
        /* new size */
        size = equil->n_element + equil->p->n[CONDENSED] + roff;
        // allocate the memory for the matrix
        matrix = (double **) malloc (sizeof(double *) * size);
        for (i = 0; i < size; i++)
          matrix[i] = (double *) malloc (sizeof(double) * (size+1));
          
        // allocate the memory for the solution vector
        sol = (double *) malloc (sizeof(double) * size);
          
        /* haven't converge yet */
        convergence_ok = false;    
      }
        
      /* reset the loop counter to compute a new equilibrium */
      k = 0;
    }
    else if (equil->verbose > 1)
    {
      printf("The solution doesn't converge\n\n");
      /* ?? */
      //remove_condensed(&size, &n_condensed, equil);
    }

    if (convergence_ok || stop)
    {
      /* when the solution have converge, we could get out of the
         main loop */
      /* if there was problem, the stop flag is set and we also get out */
      break;
    }
    
    /* suppose that it will converge the next time */
    convergence_ok = true;
    
  } /* end of main loop */
  
  if (k == iteration_max)
  {
    printf("\n");
    printf("Maximum number of %d iterations attain\n", iteration_max);
    printf("Don't thrust results.\n");
  }
  else if (stop)
  {
    printf("\n");
    printf("Problem computing equilibrium...aborted.\n");
    printf("Don't thrust results.\n");
  }
  
  printf("---------------------\n");
  printf("Results\n");
  
  print_product_composition(equil);
  printf("Temperature: \t\t % .2f K\n", equil->T);
  printf("Pressure:    \t\t % .2f atm\n", equil->P);
  
  free(sol);
  free (matrix);
  return 0;
}
  




