/* equilibrium.c  -  Responsible of the chemical equilibrium          */
/* $Id: equilibrium.c,v 1.14 2000/06/14 00:27:50 antoine Exp $ */
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <pinese@cyberwizards.com.au>                         */
/*                                                                     */
/* Licensed under the GPLv2                                            */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <time.h>

#include "print.h"
#include "equilibrium.h"
#include "load.h"
#include "libnum.h"

#include "conversion.h"
#include "compat.h"
#include "return.h"

/* Initial temperature estimate for problem with not-fixed temperature */
#define ESTIMATED_T 3800

#define MAX_PRODUCT 400
#define MAX_ELEMENT 109

/* global variable containing the information about chemical species */
propellant_t	*propellant_list;
thermo_t	    *thermo_list;


#define CONC_TOL       1.0e-8
#define LOG_CONC_TOL -18.420681
//double conc_tol   = 1.0e-8;
double conv_tol   = 0.5e-5;

int iteration_max = 400;

/* 1 for verbose, 0 for non-verbose */
int global_verbose = 0;
int global_error   = 0;

FILE * errorfile;
FILE * outputfile;

/****************************************************************
VARIABLE: Contain the molar mass of element by atomic number
          molar_mass[0] contain hydrogen and so on.
          Data come from Sargent-Welch 1996
*****************************************************************/
const float molar_mass[N_SYMB] = { 
  1.00794,   4.002602, 6.941,      9.012182, 10.811,    12.0107,
  14.00674,  15.9994,  18.9984032, 20.11797, 22.989770, 24.305, 
  26.981538, 28.0855,  30.973761,  32.066,   35.4527,   39.948, 
  39.0983,   40.078,   44.95591,   47.88,    50.9415,   51.996,
  54.938,    55.847,   58.9332,    58.6934,  63.546,    65.39,
  69.723,    72.61,    74.9216,    78.96,    79.904,    83.80,
  85.4678,   87.62,    88.9059,    91.224,   92.9064,   95.94,
  98.0,      101.07,   102.9055,   106.42,   107.868,   112.41,
  114.82,    118.71,   121.757,    127.60,   126.9045,  131.29,
  132.9054,  137.33,   138.9055,   140.12,   140.9077,  144.24,
  145.,      150.36,   151.965,    157.25,   158.9253,  162.50,
  164.9303,  167.26,   168.9342,   173.04,   174.967,   178.49,
  180.9479,  183.85,   186.207,    190.2,    192.22,    195.08,
  196.9665,  200.59,   204.383,    207.2,    208.9804,  209.,
  210.,      222.,     223.,       226.0254, 227.,      232.0381,
  231.0359,  238.029,  237.0482,   244.,     12.011,    9.01218,
  10.811,    24.305,   26.98154,   257.0,    0,         2};


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
const float  R = 8.31451; 

/* Enthalpy in the standard state */
double enthalpy_0(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);

  double val;
  int    pos = 0, i;
  
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
  val = -s->param[pos][0]*pow(T, -2) + s->param[pos][1]*pow(T, -1)*log(T)
    + s->param[pos][2] + s->param[pos][3]*T/2 + s->param[pos][4]*pow(T, 2)/3
    + s->param[pos][5]*pow(T, 3)/4 + s->param[pos][6]*pow(T, 4)/5
    + s->param[pos][7]/T;

  return val; /* dimensionless enthalpy */
}

/* Entropy in the standard state */
double entropy_0(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);
  double val;
  int    pos = 0, i;

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
    + s->param[pos][5]*pow(T, 3)/3 + s->param[pos][6]*pow(T, 4)/4 
    + s->param[pos][8];
  
  return val;
}

/* Specific heat in the standard state */
double specific_heat_0(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);
  double val;
  int    pos = 0, i;

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
}

/* 0 if out of range, 1 if ok */
int temperature_check(int sp, float T)
{
  thermo_t *s = (thermo_list + sp);

  if ((T > s->range[s->nint-1][1]) || (T < s->range[0][0]))
    return 0;

  return 1;
}


double entropy(int sp, state_t st, double ln_nj_n, float T, float P)
{
  double s;
  
  switch (st)
  {
    case GAS:
        /* The thermodynamic data are based on a standard state pressure
           of 1 bar (10^5 Pa) */
        s = entropy_0(sp, T) - ln_nj_n - log(P * ATM_TO_BAR);
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
double gibbs(int sp, state_t st, double ln_nj_n, float T, float P)
{
  double g;
  
  switch (st)
  {
    case GAS:    
        g = gibbs_0(sp, T) + ln_nj_n + log(P * ATM_TO_BAR);
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
  int i = 0, coef;
  double ans = 0;

  while ((coef = (propellant_list + molecule)->coef[i]))
	{
		ans += coef * molar_mass[(propellant_list + molecule)->elem[i]];
		i++;
	}
	return ans;
}

/* J/mol */
double heat_of_formation(int molecule)
{
  double hf = (propellant_list + molecule)->heat * 
    propellant_molar_mass(molecule);
  return hf;
}


double propellant_enthalpy(equilibrium_t *e)
{
  int i;
  double h = 0.0;
  for (i = 0; i < e->c.ncomp; i++)
  {
    h += e->c.coef[i] * heat_of_formation (e->c.molecule[i])
      / propellant_mass (e);
  }
  return h;
}

double product_enthalpy(equilibrium_t *e)
{
  int i;
  double h = 0.0;
  
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    h += e->p.coef[i][GAS] * enthalpy_0(e->p.species[i][GAS], e->T);
  }
  
  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {
    h += e->p.coef[i][CONDENSED] * enthalpy_0(e->p.species[i][CONDENSED], e->T);
  }
  return h;
}

double internal_energy(equilibrium_t *e)
{
  int i;
  double ie = 0.0;

  for (i = 0; i < e->p.n[GAS]; i++)
  {
    ie += e->p.coef[i][GAS]*enthalpy_0(e->p.species[i][GAS], e->T)
      - e->p.coef[i][GAS]/R;
  }
  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {
    ie += e->p.coef[i][CONDENSED]*enthalpy_0(e->p.species[i][CONDENSED], e->T)
      - e->p.coef[i][GAS]/R;
  }
  return ie;
}
  
double product_entropy(equilibrium_t *e)
{
  int i;
  double ent = 0.0;
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    ent += e->p.coef[i][GAS]*entropy(e->p.species[i][GAS], GAS,
                                     e->ln_nj[i] - e->ln_n,
                                     e->T, e->P);
  }
  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {
    ent += e->p.coef[i][CONDENSED]*entropy(e->p.species[i][CONDENSED],
                                           CONDENSED, 0, e->T, e->P);
  }
  return ent;
}

int thermo_search(char *str)
{
  int i;
  int last = -1;
  
  for (i = 0; i < num_thermo; i++)
  {
    if (!(STRNCASECMP(str, (thermo_list + i)->name, strlen(str))))
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
  
  for (i = 0; i < num_propellant; i++)
  {
    if (!(STRNCASECMP(str, (propellant_list + i)->name, strlen(str))))
    {
      last = i;
      printf("%-5d %s\n", i, (propellant_list + i)->name);
    }
  }
  return last;
  
}

int atomic_number(char *symbole)
{
  int i;
  int element = -1;
  
  /* find the atomic number of the element */
  for (i = 0; i < N_SYMB; i++)
  {
    if (!strcmp(symbole, symb[i]))
    {
      element = i;
      break;
    }
  }
  return element;
}


/* This fonction return the offset of the molecule in the propellant_list
   the argument is the chemical formula of the molecule */
int propellant_search_by_formula(char *str)
{
  int i = 0, j ;
  
  char  tmp[5];
  char *ptr;
  
  int   elem[6] = {0, 0, 0, 0, 0, 1};   
  int   coef[6] = {0, 0, 0, 0, 0, 0}; 

  int molecule = -1;
  
  ptr = str; /* beginning of the string */

  while ( (i < 6) && ((ptr - str) < strlen(str)) )
  {    
    if (isupper(*ptr) && islower(*(ptr+1)) && (isupper(*(ptr+2)) ||
                                               iscntrl(*(ptr+2))) )
    {
      tmp[0] = *ptr;
      tmp[1] = toupper(*(ptr+1));
      tmp[2] = '\0';
      /* find the atomic number of the element */
      elem[i] = atomic_number(tmp);
      coef[i] = 1;
      i++;   
      ptr += 2;
    }
    else if (isupper(*ptr) && (isupper(*(ptr+1)) ||
                               iscntrl(*(ptr+1))) )
    {
      tmp[0] = *ptr;
      tmp[1] = ' ';
      tmp[2] = '\0';
      elem[i] = atomic_number(tmp);
      coef[i] = 1;
      i++;
      ptr++;
    }
    else if (isupper(*ptr) && isdigit(*(ptr+1)))
    {
      tmp[0] = *ptr;
      tmp[1] = ' ';
      tmp[2] = '\0';
      elem[i] = atomic_number(tmp);
      
      j = 0;
      do
      {
        tmp[j] = *(ptr + 1 + j);
        j++;
      } while (isdigit(*(ptr + 1 + j)));

      tmp[j] = '\0';
      
      coef[i] = atoi(tmp);
      i++;
      
      ptr = ptr + j + 1;
    }
    else if (isupper(*ptr) && islower(*(ptr+1)) && isdigit(*(ptr+2)))
    {
      tmp[0] = *ptr;
      tmp[1] = toupper(*(ptr+1));
      tmp[2] = '\0';
      elem[i] = atomic_number(tmp);
      
      j = 0;
      while (isdigit(*(ptr + 2 + j)))
      {
        tmp[j] = *(ptr + 1 + j);
        j++;
      }
      tmp[j] = '\0';
      
      coef[i] = atoi(tmp);
      i++;
      
      ptr = ptr + j + 2;
    }
  }

  /*
  for (i = 0; i < 6; i++)
  {
    if (elem[i] != -1)
      printf("%s %d\n", symb[elem[i]], coef[i]);
  }
  */

  for (i = 0; i < num_propellant; i++)
  {
    for (j = 0; j < 6; j++)
    {
      /* set to the same value as the previous one if the same */
      if (!( ((propellant_list+i)->coef[j] == coef[j]) &&
             ((propellant_list+i)->elem[j] == elem[j]) ))
        break;
    }
    
  
  /* Now search in propellant list for this molecule */
/*
  for (j = 0; j < num_propellant; j++)
  {
    for (i = 0; i < 6; i++)
    {
      if ( (coef[i] != propellant_element_coef(elem[i], j)) &&
           (propellant_list + i)
        break;
    }
*/  

    if (j == 5) /* we found the molecule ! */
    {

      /* check if the inverse is true */
      molecule = i;
      break;
    }
  }
  
  return molecule;
}


double product_molar_mass(equilibrium_t *e)
{
  return (1/e->n);
}

int list_element(equilibrium_t *e)
{
  int n = 0;
  int t = 0;
  int i, j, k;

  /* reset the lement vector to -1 */
  reset_element_list(e);

  for (i = 0; i < e->c.ncomp; i++)
  {
    /* maximum of 6 different atoms in the composition */
    for (j = 0; j < 6; j++)
    {	       
      if (!( (propellant_list + e->c.molecule[i])->coef[j] == 0))
      {
        /* get the element */
        t = (propellant_list + e->c.molecule[i])->elem[j];
        
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
    fprintf(outputfile, "%d different elements in the propellant\n", n);
    /* Print those elements */
    for (i = 0; i < n; i++)
      fprintf(outputfile, "%s ", symb[e->element[i]] );
    fprintf(outputfile, "\n");
  }

  e->is_listed = 1;
  
  return n;
}

int list_product(equilibrium_t *e)
{
  int i, j, k;

  int n = 0;   /* global counter (number of species found) */
  int st;      /* temporary variable to hold the state of one specie */
  int ok = 1;

  /* reset the product to zero */
  e->p.n[GAS]       = 0;
  e->p.n[CONDENSED] = 0;
  
  for (j = 0; j < num_thermo; j++)
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

      e->p.species[ e->p.n[st] ][st] = j;
      e->p.n[st]++;
      n++;
      
      if ((e->p.n[GAS] > MAX_PRODUCT) || (e->p.n[CONDENSED] > MAX_PRODUCT))
      {
        fprintf(errorfile,
                "Error: Maximum of %d differents product reach.\n",
                MAX_PRODUCT);
        fprintf(errorfile, "       Change MAX_PRODUCT and recompile!\n");
        return ERROR;
      }
       
    }
    ok = 1;
  }
  
  /* initialize tho mol number to 0.1mol/(nb of gazeous species) */

  e->n    = e->sumn = 0.1;
  e->ln_n = log(e->n);
  
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    e->p.coef[i][GAS] = 0.1 / e->p.n[GAS];
    e->ln_nj[i]       = log( e->p.coef[i][GAS] );
  }
  
  /* initialize condensed to zero */
  for (i = 0; i < e->p.n[CONDENSED]; i++)
    e->p.coef[i][CONDENSED] = 0;

  if (e->verbose > 0)
  {
    fprintf(outputfile, "%d possible combustion product\n", n);
    fprintf(outputfile, "%d gazeous species\n", e->p.n[GAS]);
    if (e->verbose > 1)
      print_gazeous(e->p);
    fprintf(outputfile, "%d condensed species\n", e->p.n[CONDENSED]);
    if (e->verbose > 1)
      print_condensed(e->p);
  }

  return n;
}


int initialize_product(product_t *p)
{
  int i, j;
  
  for (i = 0; i < STATE_LAST; i++)
    p->n[i] = 0;
    
  /* I'm not totally sure about this typesetting - I have only allocated
     multidimensional arrays from the heap in C++ using new, never
     before in C using malloc.  However, it seems sound.
     Mark Pinese */

  /* After looking to this, it work well but we have to invert the indice
     ...[GAS][i] become ...[i][GAS]
     Antoine Lefebvre */
  if ((p->species = (int (*)[STATE_LAST]) malloc(sizeof(int) * MAX_PRODUCT *
                                                 STATE_LAST)) == NULL)
  {
	  return ERR_MALLOC;
  }
  else if ((p->coef = (double (*)[STATE_LAST]) malloc
            (sizeof(double) * MAX_PRODUCT * STATE_LAST)) == NULL)
  {
	  free(p->species);
	  return ERR_MALLOC;
  }

  p->isalloc = 1;
  
  /* initialize the list to -1 */
  for (i = 0; i < MAX_PRODUCT; i++)
  {
    for (j = 0; j < STATE_LAST; j++)
      p->species[i][j] = -1;
  }
  
  return 0;
}


int initialize_equilibrium(equilibrium_t *e)
{ 
  
  /* allocate the vector containing the delta ln(nj) */
  e->delta_ln_nj = (double *) calloc (MAX_PRODUCT, sizeof(double));
  e->ln_nj       = (double *) calloc (MAX_PRODUCT, sizeof(double));
  
  /* allocate the element vector */
  e->element = (int *) malloc (sizeof(int) * MAX_ELEMENT);
 
  /* the composition have not been set */
  e->c.ncomp = 0;
  
  e->verbose = 1;
  e->short_output = false;
  
  e->isequil = false;
  e->is_listed = 0; /* the element haven't been listed */
  
  /* initialize the product */
  return initialize_product(&(e->p));

}

int copy_equilibrium(equilibrium_t *dest, equilibrium_t *src)
{
  dest->verbose      = src->verbose;
  dest->short_output = src->short_output;
  dest->n_element    = src->n_element;
  dest->is_listed    = src->is_listed;
  dest->T            = src->T;
  dest->P            = src->P;
  dest->isequil      = src->isequil;
  dest->n            = src->n;
  dest->delta_ln_n   = src->delta_ln_n;
  dest->c.ncomp = src->c.ncomp;
  dest->p.isalloc     = src->p.isalloc;
  dest->p.n_condensed = src->p.n_condensed;
  
  memcpy(dest->delta_ln_nj, src->delta_ln_nj, MAX_PRODUCT*sizeof(double));
  memcpy(dest->element,     src->element,     MAX_ELEMENT*sizeof(int));
  memcpy(dest->c.molecule,  src->c.molecule,  MAX_COMP*sizeof(int));
  memcpy(dest->c.coef,      src->c.coef,      MAX_COMP*sizeof(double));
  memcpy(dest->p.n,         src->p.n,         STATE_LAST*sizeof(int));

  memcpy(dest->p.species,src->p.species,STATE_LAST*MAX_PRODUCT*sizeof(int));
  memcpy(dest->p.coef,   src->p.coef,   STATE_LAST*MAX_PRODUCT*sizeof(double));
  
  return 0;
}

int reset_element_list(equilibrium_t *e)
{
  int i;
  for (i = 0; i < MAX_ELEMENT; i++)
    e->element[i] = -1;
  return 0;
}

int reset_equilibrium(equilibrium_t *e)
{
  int i;

  e->n = 0.1;
  /* initialize tho mol number to 0.1mol/(nb of gazeous species) */
  for (i = 0; i < e->p.n[GAS]; i++)
    e->p.coef[i][GAS] = 0.1/e->p.n[GAS];
  
  /* initialize condensed to zero */
  for (i = 0; i < e->p.n[CONDENSED]; i++)
    e->p.coef[i][CONDENSED] = 0;

  return 0;
}

int dealloc_equilibrium(equilibrium_t *e)
{
  dealloc_product (&(e->p));
  free (e->delta_ln_nj);
  free (e->ln_nj);
  free (e->element);
  return 0;
}

int dealloc_product(product_t *p)
{
  if (!p->isalloc)
		return ERR_NOT_ALLOC;
  
  free (p->species);
  free (p->coef);
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
  e->c.molecule[ e->c.ncomp ] = sp;
  e->c.coef[ e->c.ncomp ]     = mol;
  e->c.ncomp++;
  return 0;
}

/* Mass of propellant in gram */
double propellant_mass(equilibrium_t *e)
{
  int i;
  double mass = 0.0;
  for (i = 0; i < e->c.ncomp; i++)
  {
    mass += e->c.coef[i] * propellant_molar_mass(e->c.molecule[i]);
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

/* not use for the moment */
/*
int initial_estimate(equilibrium_t *e)
{
  int i, j, mol;
  int components = e->p.n[GAS] + e->p.n[CONDENSED];
  double energy[components];

  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {
    mol = 0;
    for (j = 0; j < 5; j++)
      mol += (thermo_list + e->p.species[i][CONDENSED])->coef[j];

    energy[i] = (enthalpy_0(e->p.species[i][CONDENSED], e->T) -
                 entropy_0(e->p.species[i][CONDENSED], e->T))*R*e->T/mol;
    printf("%s \t %f \t %i\n",
           (thermo_list + e->p.species[i][CONDENSED])->name,
           energy[i], mol);
  }
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    mol = 0;
    for (j = 0; j < 5; j++)
      mol += (thermo_list + e->p.species[i][GAS])->coef[j];
    
    energy[i + e->p.n[CONDENSED]] = (enthalpy_0(e->p.species[i][GAS],
                                               e->T) -
                                      entropy_0(e->p.species[i][GAS],
                                              e->T))*R*e->T/mol;
    
    printf("%s \t %f \t %i\n",
           (thermo_list + e->p.species[i][GAS])->name,
           energy[i + e->p.n[CONDENSED]], mol);
  }
  
  return 0;
}
*/

#ifdef TRUE_ARRAY
int fill_equilibrium_matrix(double *matrix, equilibrium_t *e, problem_t P)
#else
int fill_equilibrium_matrix(double **matrix, equilibrium_t *e, problem_t P)
#endif
{

  int i, j, k;
  double tmp, mol;

  /* position of the right side dependeing on the type of problem */
  int roff = 2, size;
  
  double Mu[MAX_PRODUCT]; /* gibbs free energy for gases */
  double Ho[MAX_PRODUCT][STATE_LAST]; /* enthalpy in the standard state */
  
  if (P == TP)
    roff = 1;

  size = e->n_element + e->p.n[CONDENSED] + roff;
  
  mol = e->sumn;
  for (k = 0; k < e->p.n[GAS]; k++)
  {
    Mu[k] = gibbs(e->p.species[k][GAS], GAS, e->ln_nj[k] - e->ln_n,//log(e->n),
                  e->T, e->P);
    Ho[k][GAS] = enthalpy_0( e->p.species[k][GAS], e->T);
  }

  for (k = 0; k < e->p.n[CONDENSED]; k++)
  {
    Ho[k][CONDENSED] = enthalpy_0( e->p.species[k][CONDENSED], e->T);
  }
  
  /* fill the common part of the matrix */
  fill_matrix(matrix, e, P);
  
  /* delta ln(T) (for SP and HP only) */
  if (P != TP)
  {
    for (j = 0; j < e->n_element; j++)
    {
      tmp = 0.0;
      for (k = 0; k < e->p.n[GAS]; k++) {
        tmp += product_element_coef( e->element[j], e->p.species[k][GAS]) *
          e->p.coef[k][GAS] * Ho[k][GAS];
      }
#ifdef TRUE_ARRAY
      *(matrix + j + size * (e->n_element + e->p.n[CONDENSED] + 1)) = tmp;
#else
      matrix[j][ e->n_element + e->p.n[CONDENSED] + 1] = tmp;
#endif
    }
  }

  /* right side */
  for (j = 0; j < e->n_element; j++)
  {
    tmp = 0.0;
    
    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += product_element_coef( e->element[j], e->p.species[k][GAS]) * 
        e->p.coef[k][GAS] * Mu[k];
    
    /* b[i] */
    for (k = 0; k < STATE_LAST; k++)
      for (i = 0; i < e->p.n[k]; i++)
        tmp -= product_element_coef( e->element[j], e->p.species[i][k]) *
          e->p.coef[i][k];
    
    /* b[i]o */
    /* 04/06/2000 - division by propellant_mass(e) */
    for (i = 0; i < e->c.ncomp; i++)
      tmp += propellant_element_coef( e->element[j], e->c.molecule[i]) *
        e->c.coef[i] / propellant_mass(e);

#ifdef TRUE_ARRAY
    *(matrix + j + size * size) = tmp;
#else
    matrix[j][size] = tmp; 
#endif
  }

  /* delta ln(T) */
  if (P != TP)
  {
    for (j = 0; j < e->p.n[CONDENSED]; j++) /* row */
#ifdef TRUE_ARRAY
      *(matrix + j + e->n_element + size *
        (e->n_element + e->p.n[CONDENSED] + 1)) = Ho[j][CONDENSED];
#else
    matrix[ j + e->n_element ][ e->n_element + e->p.n[CONDENSED] + 1] =
      Ho[j][CONDENSED];
#endif
    
  }
  
  /* right side */
  for (j = 0; j < e->p.n[CONDENSED]; j++) /* row */
  {
#ifdef TRUE_ARRAY
    *(matrix + j + e->n_element + size * size) =
      gibbs( e->p.species[j][CONDENSED], CONDENSED, 0, e->T, e->P);
#else
    matrix[ j + e->n_element ][size] =
      gibbs( e->p.species[j][CONDENSED], CONDENSED, 0, e->T, e->P);
#endif
    
  }

  /* delta ln(n) */
#ifdef TRUE_ARRAY
  *(matrix + e->n_element + e->p.n[CONDENSED] + size *
    (e->n_element + e->p.n[CONDENSED])) = mol - e->n;
#else
  matrix[e->n_element + e->p.n[CONDENSED]][e->n_element + e->p.n[CONDENSED]]
    =  mol - e->n;
#endif
  
  /* delta ln(T) */
  if (P != TP)
  {
    tmp = 0.0;
    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += e->p.coef[k][GAS] * Ho[k][GAS];
#ifdef TRUE_ARRAY
    *(matrix + e->n_element + e->p.n[CONDENSED] + size *
      (e->n_element + e->p.n[CONDENSED] + 1)) = tmp;
#else
    matrix[e->n_element + e->p.n[CONDENSED]][e->n_element
                                             + e->p.n[CONDENSED]+ 1] = tmp;
#endif
    
  }
  
  /* right side */
  tmp = 0.0;
  for (k = 0; k < e->p.n[GAS]; k++)
  {
    tmp += e->p.coef[k][GAS] * Mu[k];
  }

#ifdef TRUE_ARRAY
  *(matrix + e->n_element + e->p.n[CONDENSED] + size * size) =
    e->n - mol + tmp;
#else
  matrix[e->n_element + e->p.n[CONDENSED] ][size] = e->n - mol + tmp;
#endif
  
  
  /* for enthalpy/pressure problem */
  if (P == HP)
  {
    /* part with lagrangian multipliers */
    for (i = 0; i < e->n_element; i++) /* each column */
    {   
      tmp = 0.0;
      for (k = 0; k < e->p.n[GAS]; k++)
        tmp += product_element_coef( e->element[i], e->p.species[k][GAS] ) * 
          e->p.coef[k][GAS] * Ho[k][GAS];

#ifdef TRUE_ARRAY
      *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size * i) = tmp;
#else
      matrix[ e->n_element + e->p.n[CONDENSED] + 1 ][i] = tmp;      
#endif
    }

    /* Delta n */
    for (i = 0; i < e->p.n[CONDENSED]; i++)
#ifdef TRUE_ARRAY
      *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size *
        ( i + e->n_element)) = Ho[i][CONDENSED];
#else
      matrix[ e->n_element + e->p.n[CONDENSED] + 1 ][i + e->n_element] = 
        Ho[i][CONDENSED];
#endif

    /* Delta ln(n) */
    tmp = 0.0;
    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += e->p.coef[k][GAS] * Ho[k][GAS];

#ifdef TRUE_ARRAY
    *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size *
      (e->n_element + e->p.n[CONDENSED])) = tmp;
#else
    matrix[e->n_element + e->p.n[CONDENSED] + 1][e->n_element + 
                                                 e->p.n[CONDENSED] ] = tmp;
#endif

    /* Delta ln(T) */
    tmp = 0.0;
    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += e->p.coef[k][GAS]*specific_heat_0( e->p.species[k][GAS], e->T );

    for (k = 0; k < e->p.n[CONDENSED]; k++)
      tmp += e->p.coef[k][CONDENSED]*
        specific_heat_0( e->p.species[k][CONDENSED], e->T);

    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += e->p.coef[k][GAS] * Ho[k][GAS] * Ho[k][GAS];

#ifdef TRUE_ARRAY
    *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size *
      (e->n_element + e->p.n[CONDENSED] + 1)) = tmp;
#else
    matrix[e->n_element + e->p.n[CONDENSED] + 1][e->n_element +
                                                 e->p.n[CONDENSED] + 1] = tmp;
#endif
    
    
    /* right side */
    tmp = 0.0;
    tmp = propellant_enthalpy(e)/(R*e->T) - product_enthalpy(e);
    
    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += e->p.coef[k][GAS] * Ho[k][GAS] * Mu[k];

#ifdef TRUE_ARRAY
    *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size * size) = tmp;
#else
    matrix[e->n_element + e->p.n[CONDENSED] + 1][size] = tmp;
#endif
    
  } /* for entropy/pressure problem */
  else if (P == SP)
  {
    /* part with lagrangian multipliers */
    for (i = 0; i < e->n_element; i++) /* each column */
    {   
      tmp = 0.0;
      for (k = 0; k < e->p.n[GAS]; k++)
        tmp += product_element_coef( e->element[i], e->p.species[k][GAS] ) * 
          e->p.coef[k][GAS] *
          entropy(e->p.species[k][GAS], GAS,
                  e->ln_nj[k] - e->ln_n,//log(e->n),
                  e->T, e->P);
      
#ifdef TRUE_ARRAY
      *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size * i) = tmp;
#else
      matrix[ e->n_element + e->p.n[CONDENSED] + 1][i] = tmp;      
#endif 
    }
    
    /* Delta n */
    for (i = 0; i < e->p.n[CONDENSED]; i++)
#ifdef TRUE_ARRAY
      *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size *
        (i + e->n_element)) = entropy_0( e->p.species[i][CONDENSED], e->T);
#else
      matrix[ e->n_element + e->p.n[CONDENSED] + 1 ][i + e->n_element] = 
        entropy_0( e->p.species[i][CONDENSED], e->T); /* ok for condensed */
#endif
    

    /* Delta ln(n) */
    tmp = 0.0;
    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += e->p.coef[k][GAS] * entropy(e->p.species[k][GAS], GAS,
                                         e->ln_nj[k] - e->ln_n,//log(e->n),
                                         e->T, e->P);

#ifdef TRUE_ARRAY
    *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size *
      (e->n_element + e->p.n[CONDENSED] )) = tmp;
#else
    matrix[e->n_element + e->p.n[CONDENSED] + 1][e->n_element + 
                                                 e->p.n[CONDENSED] ] = tmp;
#endif
    
    tmp = 0.0;
    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += e->p.coef[k][GAS]*specific_heat_0( e->p.species[k][GAS], e->T );

    for (k = 0; k < e->p.n[CONDENSED]; k++)
      tmp += e->p.coef[k][CONDENSED]*
        specific_heat_0( e->p.species[k][CONDENSED], e->T);

    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += e->p.coef[k][GAS]*Ho[k][GAS]*
        entropy(e->p.species[k][GAS], GAS,
                e->ln_nj[k] - e->ln_n,//log(e->n),
                e->T, e->P);
    
#ifdef TRUE_ARRAY
    *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size *
      (e->n_element + e->p.n[CONDENSED] + 1)) = tmp;
#else
    matrix[e->n_element + e->p.n[CONDENSED] + 1][e->n_element +
                                                 e->p.n[CONDENSED] + 1] = tmp;
#endif
    
    
    /* entropy of reactant */
    tmp = e->entropy; /* assign entropy */
    tmp -= product_entropy(e);
    tmp += e->n;

    for (k = 0; k < e->p.n[GAS]; k++)
      tmp -= e->p.coef[k][GAS];

    for (k = 0; k < e->p.n[GAS]; k++)
      tmp += e->p.coef[k][GAS]
        * Mu[k]
        * entropy(e->p.species[k][GAS], GAS,
                  e->ln_nj[k] - e->ln_n,//log(e->n),
                  e->T, e->P);

#ifdef TRUE_ARRAY
    *(matrix + e->n_element + e->p.n[CONDENSED] + 1 + size * size) = tmp;
#else
    matrix[e->n_element + e->p.n[CONDENSED] + 1][size] = tmp;
    
#endif
    
  }

  return 0;
}

/* This part of the matrix is the same for equilibrium and derivative */
#ifdef TRUE_ARRAY
int fill_matrix(double *matrix, equilibrium_t *e, problem_t P)
#else
int fill_matrix(double **matrix, equilibrium_t *e, problem_t P)
#endif
{

  int i, j, k, size;
  double tmp;

  int roff = 2;

  if (P == TP)
    roff = 1;
  
  size = e->n_element + e->p.n[CONDENSED] + roff;
  
  /* fill the matrix (part with the Lagrange multipliers) */
  for (i = 0; i < e->n_element; i++)  /* each column */
  {
    for (j = 0; j < e->n_element; j++) /* each row */
    {
      tmp = 0.0;
      for (k = 0; k < e->p.n[GAS]; k++)
      {
        tmp += product_element_coef( e->element[j], e->p.species[k][GAS]) * 
          product_element_coef( e->element[i], e->p.species[k][GAS] ) * 
          e->p.coef[k][GAS]; 
      }
#ifdef TRUE_ARRAY
      *(matrix + j + size * i) = tmp;
#else
      matrix[j][i] = tmp;
#endif 
    }
  }
  
  /* Delta n */
  for (i = 0; i < e->p.n[CONDENSED]; i++) /* column */
  {
    for (j = 0; j < e->n_element; j++) /* row */
    {
#ifdef TRUE_ARRAY
      *(matrix + j + size * (i + e->n_element)) =
        product_element_coef(e->element[j], e->p.species[i][CONDENSED]);
#else
      matrix[j][i + e->n_element ] =
        product_element_coef(e->element[j], e->p.species[i][CONDENSED]);
#endif 
    }
  } 

  /* delta ln(n) */
  for (j = 0; j < e->n_element; j++)
  {
    tmp = 0.0;
    for (k = 0; k < e->p.n[GAS]; k++) {
      tmp += product_element_coef( e->element[j], e->p.species[k][GAS]) * 
        e->p.coef[k][GAS];
    }
#ifdef TRUE_ARRAY
    *(matrix + j + size * (e->n_element + e->p.n[CONDENSED])) = tmp;
#else
    matrix[j][ e->n_element + e->p.n[CONDENSED] ] = tmp;
#endif
  }
   
  /* second row */
  for (i = 0; i < e->n_element; i++) /* column */
  {
    for (j = 0; j < e->p.n[CONDENSED]; j++) /* row */
    {
      /* copy the symetric part of the matrix */
#ifdef TRUE_ARRAY
      *(matrix + j + e->n_element + size * i) =
        *(matrix + i + size * (j + e->n_element));
#else
      matrix[j + e->n_element ][i] = matrix[i][j + e->n_element ];
#endif 
    }
  }
  
  /* set to zero */
  for (i = 0; i < e->p.n[CONDENSED]+1; i++) /* column */
  {
    for (j = 0; j < e->p.n[CONDENSED]; j++) /* row */
    {
#ifdef TRUE_ARRAY
      *(matrix + j + e->n_element + size * (i + e->n_element)) = 0.0;
#else
      matrix[j + e->n_element ][i + e->n_element] = 0.0;
#endif 
    }
  }
  
  /* third row */
  for (i = 0; i < e->n_element; i++) /* each column */
  {   
    /* copy the symetric part of the matrix */
#ifdef TRUE_ARRAY
    *(matrix + e->n_element + e->p.n[CONDENSED] + size * i) =
      *(matrix + i + size * (e->n_element + e->p.n[CONDENSED]));
#else
    matrix[ e->n_element + e->p.n[CONDENSED] ][i] =
      matrix[i][ e->n_element + e->p.n[CONDENSED] ];
#endif
  }

  /* set to zero */
  for (i = 0; i < e->p.n[CONDENSED]; i++) /* column */
  {
#ifdef TRUE_ARRAY
    *(matrix + e->n_element + e->p.n[CONDENSED] + size *
      (i + e->n_element)) = 0.0;
#else
      matrix[e->n_element + e->p.n[CONDENSED] ][i + e->n_element] = 0.0;
#endif 
  }
  
  return 0;
}



/* may be optimize!!!!!!!! */
int remove_condensed(int *size, int *n, equilibrium_t *e)
{

  int i, j, k, pos;
  int r = 0; /* something have been replace, 0=false, 1=true */

  int ok = 1;
  
  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {

    /* if a condensed have negative coefficient, we should remove it */
    if (e->p.coef[i][CONDENSED] <= 0.0)
    {
      if (e->verbose > 1)
      {
        fprintf(outputfile, "%s should be remove\n", 
                (thermo_list + e->p.species[i][CONDENSED])->name );
      }
      /* remove from the list ( put it at the end for later use )*/
      pos = e->p.species[i][CONDENSED];
      
      for (j = i; j < e->p.n[CONDENSED] - 1; j++)
      {
        e->p.species[j][CONDENSED] = 
          e->p.species[j + 1][CONDENSED];
      }
      e->p.species[ e->p.n[CONDENSED] - 1 ][CONDENSED] = pos;
        
      e->p.n[CONDENSED] = e->p.n[CONDENSED] - 1;
      (*size)--; /* reduce the size of the matrix */
      r = 1;
    }
    else if ( !(temperature_check(e->p.species[i][CONDENSED], e->T)) )
    {
      /* if the condensed species is present outside of the temperature
         range at which it could exist, we should replace it by an other
         phase */

      /* Find the new molecule */
      for (j = e->p.n[CONDENSED]; j < (*n); j++)
      {
        /* if this is the same molecule and temperature_check is true,
           than replace the molecule */

        for (k = 0; k < 5; k++)
        {
          if (!( ((thermo_list + e->p.species[i][CONDENSED])->coef[k] ==
                  (thermo_list + e->p.species[j][CONDENSED])->coef[k] ) &&
                 ((thermo_list + e->p.species[i][CONDENSED])->elem[k] ==
                  (thermo_list + e->p.species[j][CONDENSED])->elem[k] ) &&
                 temperature_check(e->p.species[j][CONDENSED], e->T) ))
          {
            ok = 0;
          }
        }

        /* replace the molecule */
        if (ok)
        {
          if (e->verbose > 1)
          {
            fprintf(outputfile, "%s should be replace by %s\n",
                    (thermo_list + e->p.species[i][CONDENSED])->name,
                    (thermo_list + e->p.species[j][CONDENSED])->name);
          }

          pos = e->p.species[i][CONDENSED];
          e->p.species[i][CONDENSED] = e->p.species[j][CONDENSED];
          e->p.species[j][CONDENSED] = pos;

          r = 1; /* A species have been replace */
          
          /* we do not need to continue searching so we break */
          break;
        }

        ok = 1;
      }
    }
  } /* for each condensed */
  
  /* 0 if none remove */
  return r;
}

int include_condensed(int *size, int *n, equilibrium_t *e, 
                      double *sol)
{
  double tmp;
  double temp;
  int    i, j, k;
  int    pos;

  tmp = 0.0;
  j   = -1;

  /* We include a condensed if it minimize the gibbs free energy and
     if it could exist at the chamber temperature */
  for (i = e->p.n[CONDENSED] ; i < (*n); i++)
  {
    if (temperature_check(e->p.species[i][CONDENSED], e->T))
    {
      temp = 0.0;
      for (k = 0; k < e->n_element; k++)
        temp += sol[k]*product_element_coef(e->element[k], 
                                            e->p.species[i][CONDENSED]);
      
      if ( gibbs_0( e->p.species[i][CONDENSED], e->T) - temp < tmp )
      {
        tmp = gibbs_0( e->p.species[i][CONDENSED], e->T) - temp;
        j = i; 
      }
    }
  }

  /* In the case we found a species that minimize the gibbs energy,
     we should include it */
  if (!(j == -1))
  {
    if (e->verbose > 1)
    {
      fprintf(outputfile, "%s should be include\n", 
              (thermo_list + e->p.species[j][CONDENSED])->name );
    }
    /* to include the species, exchange the value */
    pos = e->p.species[ e->p.n[CONDENSED] ][CONDENSED];
    e->p.species[ e->p.n[CONDENSED] ][CONDENSED] = 
      e->p.species[j][CONDENSED];
    e->p.species[j][CONDENSED] = pos;
    
    e->p.n[CONDENSED]++;
    
    return 1;
  }
  return 0;
}


int new_approximation(equilibrium_t *e, double *sol, problem_t P)
{
  int i, j;

  /* control factor */
  double lambda1, lambda2, lambda;
  //double delta_ln_T;
  
  double temp;
    
  /* compute the values of delta ln(nj) */
  e->delta_ln_n = sol[ e->n_element + e->p.n[CONDENSED] ];

  if  (P != TP)
    e->delta_ln_T = sol[e->n_element + e->p.n[CONDENSED] + 1];
  else
    e->delta_ln_T = 0.0;

  
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    temp = 0.0;
    for (j = 0; j < e->n_element; j++)
    {
      temp +=
        product_element_coef(e->element[j], e->p.species[i][GAS]) * sol[j];
    }
    
    e->delta_ln_nj[i] =
      - gibbs(e->p.species[i][GAS], GAS, e->ln_nj[i] - e->ln_n,//log(e->n),
              e->T, e->P)
      + temp
      + e->delta_ln_n
      + enthalpy_0(e->p.species[i][GAS], e->T)*e->delta_ln_T;     
  }
  

  lambda2 = 1.0;
  lambda1 = __max(fabs(e->delta_ln_T), fabs(e->delta_ln_n));
  lambda1 = 5 * lambda1;
  
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    if (e->delta_ln_nj[i] > 0.0)
    {
      if (e->ln_nj[i] - e->ln_n /*log(e->n)*/ <= LOG_CONC_TOL)//log(conc_tol)) 
      {

        lambda2 = __min(lambda2,
                        fabs( ((- e->ln_nj[i] + e->ln_n /*log(e->n)*/
                                - 9.2103404)
                               /(e->delta_ln_nj[i] - e->delta_ln_n))) );
//        lambda2 = __min(lambda2,
//                        fabs( ((-log(e->p.coef[i][GAS]/e->n) - 9.2103404)
//                               /(e->delta_ln_nj[i] - e->delta_ln_n))) );
      }
      else if (e->delta_ln_nj[i] > lambda1) 
      {
        lambda1 = e->delta_ln_nj[i];
      }
      
    }
  }

  lambda1 = 2.0 / lambda1;
 
  lambda = _min(1.0, lambda1, lambda2);
  
  if (e->verbose > 2)
  {
    fprintf(outputfile, "lambda = %.10f, lambda1 = %.10f, lambda2 = %.10f\n",
            lambda, lambda1, lambda2);
    fprintf(outputfile, " \t  nj \t\t  ln_nj_n \t Delta ln(nj)\n");

    for (i = 0; i < e->p.n[GAS]; i++)
    {
      if (1)//!(e->p.coef[i][GAS] == 0))
        fprintf(outputfile, "%s \t % .4e \t % .4e \t % .4e\n", 
                (thermo_list + e->p.species[i][GAS])->name, 
                e->p.coef[i][GAS],
                e->ln_nj[i],
                e->delta_ln_nj[i]);
    }
  }

  e->sumn = 0.0;
  
  /* compute the new value for nj (gazeous) and ln_nj */
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    e->ln_nj[i] = e->ln_nj[i] + lambda * e->delta_ln_nj[i];

    if ( (e->ln_nj[i] - e->ln_n /*log(e->n)*/) <= LOG_CONC_TOL)//log(conc_tol))
    {
      e->p.coef[i][GAS] = 0.0;
    }
    else
    {
      e->p.coef[i][GAS] = exp(e->ln_nj[i]);
      e->sumn += e->p.coef[i][GAS];
    }
    
  }
  
  /* compute the new value for nj (condensed) */
  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {
    e->p.coef[i][CONDENSED] = e->p.coef[i][CONDENSED] + 
      lambda*sol[e->n_element + i];     
  }

  if (e->verbose > 2)
  {
    for (i = 0; i < e->p.n[CONDENSED]; i++)
    {
      fprintf(outputfile, "%s: \t %f\n", 
              (thermo_list + e->p.species[i][CONDENSED])->name, 
              e->p.coef[i][CONDENSED]);
    }
  }
    
  /* new value of T */
  if (P != TP)
    e->T = exp( log(e->T) + lambda * e->delta_ln_T);
      
  if (e->verbose > 2)
    fprintf(outputfile, "Temperature: %f\n", e->T);
      
  /* new value of n */
  e->ln_n = e->ln_n + lambda * e->delta_ln_n;
  e->n = exp(e->ln_n);
  
  //e->n = exp (log(e->n) + lambda * e->delta_ln_n);
  
  return SUCCESS;
}
    
bool convergence(equilibrium_t *e, double *sol)
{
  int i;
  double mol;

  /* for convergence test, mol is the sum of all mol
     even condensed */
  
  mol = e->sumn;
      
  /* check for convergence */ 
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    if (!(e->p.coef[i][GAS]*fabs(e->delta_ln_nj[i])/mol <= conv_tol))
      return false; /* haven't converge yet */
  }
      
  for ( i = 0; i < e->p.n[CONDENSED]; i++ )
  {
    /* test for the condensed phase */
    if (!(sol[e->n_element+1]/mol <= conv_tol))
      return false; /* haven't converge yet */
  }
      
  if (!(e->n*fabs(e->delta_ln_n)/mol <= conv_tol))
    return false; /* haven't converge yet */

  if (!(fabs(e->delta_ln_T) <= 1.0e-4))
    return false;
  
  return true;
}

int equilibrium(equilibrium_t *equil, problem_t P)
{ 
  int       i, k;
  int       size;     /* size of the matrix */
#ifdef TRUE_ARRAY
  double *matrix;
#else
  double ** matrix;  
#endif
  double  * sol;

  bool      convergence_ok;
  bool      stop           = false;
  bool      gas_reinserted = false;
  bool      solution_ok    = false;
  
  /* position of the right side of the matrix dependeing on the
     type of problem */
  int roff = 2;

  if (P == TP)
    roff = 1;

  /* initial temperature for assign enthalpy, entropy/pressure */
  if ( P != TP)
    equil->T = ESTIMATED_T;

  
  if (!(equil->is_listed)) /* if the element and the product haven't
                              been listed */
  {
    list_element(equil);

    if (list_product(equil) == ERROR)
    {
      return ERROR;
    }
    equil->p.n_condensed = equil->p.n[CONDENSED];
  }
  
  /* First determine an initial estimate of the composition
     to accelerate the convergence */
  /* initial_estimate(equil); */
  
  /* For the first equilibrium, we do not consider the condensed
     species. */
  if (!(equil->isequil))
  {
    equil->p.n[CONDENSED] = 0;
    equil->n = 0.1; /* initial estimate of the mol number */
  }
  
  /* the size of the coefficient matrix */
  size = equil->n_element + equil->p.n[CONDENSED] + roff;
  
  /* allocate the memory for the matrix */
#ifdef TRUE_ARRAY
  matrix = (double *) malloc (size*(size+1)*sizeof(double));
#else
  matrix = (double **) malloc (sizeof(double *) * size);
  for (i = 0; i < size; i++)
    matrix[i] = (double *) malloc (sizeof(double) * (size+1));
#endif
  
  /* allocate the memory for the solution vector */
  sol = (double *) calloc (size, sizeof(double));

  /* main loop */
  for (k = 0; k < iteration_max; k++)
  {
    /* Initially we haven't a good solution */
    solution_ok = false;

    while (!solution_ok)
    {      
      fill_equilibrium_matrix(matrix, equil, P);
      
      if (equil->verbose > 2)
      {
        fprintf(outputfile, "Iteration %d\n", k+1);
        print_matrix(matrix, size);
      }
      if ( lu(matrix, sol, size) == -1) /* solve the matrix */
      {
        /* the matrix have no unique solution */
        fprintf(outputfile,
                "The matrix is singular, removing excess condensed.\n");
          
        /* Try removing excess condensed */
        if (!remove_condensed(&size, &(equil->p.n_condensed), equil))
        {
          if (gas_reinserted)
          {
            fprintf(errorfile, "ERROR: No convergence, don't trust results\n");
            /* finish the main loop */
            stop = true;
            break;
          }
          fprintf(errorfile, "None remove. Try reinserting remove gaz\n");
          for (i = 0; i < equil->p.n[GAS]; i++)
          {
            /* It happen that some species were eliminated in the
               process even if they should be prsent in the equilibrium.
               In such case, we have to reinsert them */
            if (equil->p.coef[i][GAS] == 0.0)
              equil->p.coef[i][GAS] = 1e-6;
          }
          gas_reinserted = true;
        }
        else
          gas_reinserted = false;
          
        /* Restart the loop counter to zero for a new loop */
        k = 0;
      }
      else /* There is a solution */
      {
        solution_ok = true;
      }
    }
      
    if (equil->verbose > 2)
      print_vec(sol, size);    /* print the solution vector */

    /* compute the new approximation */
    new_approximation(equil, sol, P);

    convergence_ok = false;

    /* verify the convergence */
    if (convergence(equil, sol))
    {
      convergence_ok = true;

      if (equil->verbose > 0)
      {
        fprintf(outputfile, "Convergence: %-2d iteration, %f deg K\n",
                k+1, equil->T);
        //fprintf(outputfile, "T = %f\n", equil->T);
      }
      gas_reinserted = false;
      
      /* find if a new condensed species should be include or remove */
      if (remove_condensed(&size, &(equil->p.n_condensed), equil) ||
          include_condensed(&size, &(equil->p.n_condensed), equil, sol))
      {

#ifndef TRUE_ARRAY
        for (i = 0; i < size; i++)
          free(matrix[i]);
#endif
        free(matrix);

        
        free(sol);
        /* new size */
        size = equil->n_element + equil->p.n[CONDENSED] + roff;

        /* allocate the memory for the matrix */
#ifdef TRUE_ARRAY
        matrix = (double *) malloc(size*(size+1)*sizeof(double));
#else
        matrix = (double **) malloc (sizeof(double *) * size);
        for (i = 0; i < size; i++)
          matrix[i] = (double *) malloc (sizeof(double) * (size+1));
#endif
                
        /* allocate the memory for the solution vector */
        sol = (double *) malloc (sizeof(double) * size);
          
        /* haven't converge yet */
        convergence_ok = false;    
      }
        
      /* reset the loop counter to compute a new equilibrium */
      k = 0;
    }
    else if (equil->verbose > 2)
    {
      fprintf(outputfile, "The solution doesn't converge\n\n");
      /* ?? */
      /*remove_condensed(&size, &n_condensed, equil); */
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

  free (sol);

#ifndef TRUE_ARRAY
  for (i = 0; i < size; i++)
    free (matrix[i]);
#endif
  free (matrix);
  
  if (k == iteration_max)
  {
    fprintf(outputfile, "\n");
    fprintf(outputfile, "Maximum number of %d iterations attain\n",
            iteration_max);
    fprintf(outputfile, "Don't thrust results.\n"); 
    return ERROR;
  }
  else if (stop)
  {
    fprintf(outputfile, "\n");
    fprintf(outputfile, "Problem computing equilibrium...aborted.\n");
    fprintf(outputfile, "Don't thrust results.\n");
    return ERROR;
  }

  equil->isequil = true;
  
  return SUCCESS;
}







