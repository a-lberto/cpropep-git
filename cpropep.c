/* cpropep.c  -  Calculation of Complex Chemical Equilibrium           */
/* Copyright (C) 2000                                                  */
/* Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                      */
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

#include "cpropep.h"
#include "load.h"
#include "libnum.h"


/* global variable containing the information about chemical species */
propellant_t	*propellant_list;
thermo_t	*thermo_list;


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
const float  R = 8.314; 



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

      return (val * R * T); /* to convert from dimensionless to J/mol */
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
      
      return val * R; /* to convert from dimensionless to J/(mol K) */
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
      
      return val * R; /* to convert from dimensionless to J/(mol K) */
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

/* enthalpy variation between 298.15 and T */
double delta_enthalpy(int sp, float T)
{
  /* delta henthalpy in J/mol */
  return enthalpy(sp, T) - (thermo_list + sp)->heat;
}


/* J/mol T is in K, P is in atm */
double gibbs(int sp, state_t st, double nj, double n, float T, float P)
{

  double g;

  switch (st)
  {
  case GAS:    
    g = ( enthalpy(sp, T) - T * entropy(sp, T) )
      + R*T*log(nj/n) + R*T*log(P);
    break;
  case CONDENSED:
    g = ( enthalpy(sp, T) - T * entropy(sp, T) );
    
  default:
    g = 0;
  }
  
  return g/(R*T);
}


/* kJ/mol */
//double potential(int sp, float T)
//{
  /* u = Hf - TS */
//  return enthalpy(sp, T)/1000  - T *  entropy(sp, T)/1000;
//}


/* give the molar mass of a propellant or pruduct in g/mol */
/* for each element int the molecule, coefficient * molar mass of */
/* the element */
double propellant_molar_mass(int molecule)
{     
  int i;
  double ans = 0;
  for (i = 0; i < 6; i++)
    ans += (propellant_list + molecule)->coef[i] * molar_mass[(propellant_list + molecule)->elem[i]];
  return ans;
}

/* give the heat of formation of a propellant in kJ/mol */
double heat_of_formation(int molecule)
{
  return (propellant_list + molecule)->heat * 4.1868 * propellant_molar_mass(molecule) / 1000;
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
  {
    printf("%d   %s\n", i, (thermo_list + i)->name);
  }
  return 0;
}

int print_propellant_list(void)
{
  int i;
  for (i = 0; i < MAX_PROPELLANT; i++)
  {
    printf("%d   %s\n", i, (propellant_list + i)->name);
  }
  return 0;
}


int list_element(composition_t comp, int *list)
{
  int e = 0;
  int t = 0;
  int i, j, k;
  
  for (i = 0; i < comp.ncomp; i++)
  {
    for (j = 0; j < 6; j++)
    {	       
      if (!(propellant_list[ comp.molecule[i] ].coef[j] == 0))
      {
	t = propellant_list[ comp.molecule[i] ].elem[j];
	for (k = 0; k <= e; k++)
	{
	  if (list[k] == t)
	    break;
	  if (k == e)
	  {
	    list[e] = t;
	    e++;
	    break;
	  }
	}
      }
    }
  }
  return e;
}

int list_product(int n_element, int *element_list, product_t *p, float T)
{
  int i, j, k;

  int e = 0;   // global counter (number of species found)
  int st;      // temporary variable to hold the state of one specie
  int ok = 1;
  
  for (j = 0; j < MAX_THERMO; j++)
  {
    for (k = 0; k < 5; k++)
    {
      if (!(thermo_list[j].coef[k] == 0))
      {
	for (i = 0; i < n_element; i++)
	{
	  if (element_list[i] == thermo_list[j].elem[k])
	    break;
	  else if (i == n_element -1)
	    ok = 0;
	}    
	if (!ok)
	  break;
      }
    }
    if (ok) /* add to the list */
    {
      st = thermo_list[j].state;

      if (temperature_check(j, T) )  // if the molecule could exist at that temperature
      {
	p->species[st][ p->n[st] ] = j;
	p->n[st]++;
	e++;
      }

    }
    ok = 1;
  }

  // reallocate the momory to the minimum size
  for (i = 0; i < STATE_LAST; i++)
  {
    if (( p->species[i] = (int *) realloc( p->species[i], 
					   sizeof(int) * p->n[i])) == NULL)
      printf("Reallocation of memory failed\n");
    
    if (( p->coef[i] = (double *) realloc ( p->coef[i], 
					   sizeof(double) * p->n[i])) == NULL)
      printf("Reallocation of memory failed\n");
  }

  
  for (i = 0; i < p->n[GAS]; i++)
    p->coef[GAS][i] = 0.1/p->n[GAS];
 
  for (i = 0; i < p->n[CONDENSED]; i++)
    p->coef[CONDENSED][i] = 0;

  return e;
}


// find a first approximation of the composition
int approx(int n_element, int *element_list, product_t *p, composition_t *c)
{

  // doesn't work for the moment.......
  int i, j;

  double mol_i = 0.0;
  double mol_o = 0.0;

  for (i = 0; i < c->ncomp; i++)
    mol_i += c->coef[i];

  for (i = 0; i < STATE_LAST; i++)
    for (j = 0; j < p->n[i]; j++)
      mol_o += p->coef[i][j];

  j = (int)GAS;
  for (i = 0; i < (p->n[GAS] + p->n[CONDENSED] - n_element); i++)
    {
    }
  return 0;

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


/* use the theory explain in
Theorical Evaluation of chemical propellant
by Roger Lawrence Wilkins
Prentice-Hall International series in space technology */
int fill_matrix(double **matrix, int n_elements, int *element, product_t *p, 
		composition_t *b)
{

  int i, j, k;
  double tmp, mol;

// fill the matrix (part with the Lagrange multipliers)
  for (i = 0; i < n_elements; i++) // each column
  {
    for (j = 0; j < n_elements; j++) // each row
    {
      tmp = 0.0;
      for (k = 0; k < p->n[GAS]; k++)
  	tmp += product_element_coef( element[j], p->species[GAS][k]) * 
	  product_element_coef( element[i], p->species[GAS][k] ) * 
	  p->coef[GAS][k];  
      matrix[j][i] = tmp;
    }
  }
  
  // X1c
  for (i = 0; i < p->n[CONDENSED]; i++) // column
  {
    for (j = 0; j < n_elements; j++) // row
    {
      matrix[j][i + n_elements ] = product_element_coef(element[j], 
						      p->species[CONDENSED][i]);
    }
  } 

  // delta ln(n)
  for (j = 0; j < n_elements; j++)
  {
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
  	tmp += product_element_coef( element[j], p->species[GAS][k]) * 
	  p->coef[GAS][k];
    matrix[j][ n_elements + p->n[CONDENSED] ] = tmp;
  }



  mol = 0.0;
  for (k = 0; k < p->n[GAS]; k++)
    mol += p->coef[GAS][k];

  // right side
  for (j = 0; j < n_elements; j++)
  {
    tmp = 0.0;
    for (k = 0; k < p->n[GAS]; k++)
  	tmp += product_element_coef( element[j], p->species[GAS][k]) * 
	  p->coef[GAS][k] * gibbs( p->species[GAS][k], GAS, p->coef[GAS][k],
				  mol, 1200, 136);


    // b[i]
    for (k = 0; k < STATE_LAST; k++)
      for (i = 0; i < p->n[k]; i++)
	tmp -= product_element_coef( element[j], p->species[k][i]) *
	  p->coef[k][i];
    
    // b[i]o
    for (i = 0; i < b->ncomp; i++)
      tmp += propellant_element_coef( element[j], b->molecule[i]) *
	b->coef[i];


    matrix[j][ n_elements + p->n[CONDENSED] + 1 ] = tmp;  
  }


  //second row
  for (i = 0; i < n_elements; i++) // column
  {
    for (j = 0; j < p->n[CONDENSED]; j++) // row
    {
      //      printf("%d\n",j);
      matrix[j + n_elements ][i] = product_element_coef(element[i],
					       p->species[CONDENSED][j]);
    }
  }

  
  // right side
  for (j = 0; j < p->n[CONDENSED]; j++) // row
  {
    matrix[ j + n_elements ][ n_elements + p->n[CONDENSED] + 1] = 
      gibbs( p->species[GAS][j], GAS, p->coef[GAS][j],
	     mol, 1200, 136);
  
  }


  // third big row
  for (i = 0; i < n_elements; i++) // each column
  {   
      tmp = 0.0;
      for (k = 0; k < p->n[GAS]; k++)
  	tmp += product_element_coef( element[i], p->species[GAS][k] ) * 
	  p->coef[GAS][k];

      matrix[ n_elements + p->n[CONDENSED] ][i] = tmp;      
  }

  // delta ln(n)
  tmp = 0.0;
  for (k = 0; k < p->n[GAS]; k++)
    tmp += p->coef[GAS][k];
  matrix[ n_elements + p->n[CONDENSED] ][ n_elements + p->n[CONDENSED] ] = 
    tmp - mol;


  tmp = 0.0;
  for (k = 0; k < p->n[GAS]; k++)
  {
    tmp += p->coef[GAS][k]*gibbs( p->species[GAS][j], GAS, p->coef[GAS][j],
				 mol, 1200, 136) +
      p->coef[GAS][k];
  }

  matrix[n_elements + p->n[CONDENSED] ][ n_elements + p->n[CONDENSED] + 1] = 
    mol - tmp;


  return 0;
}


/* to be wrtite */
int new_approx(double *matrixsol, int n_element, int *element, product_t p ) 
{ 
  double mol_x;
  double mol_y;

  return 0;
}


int main(int argc, char *argv[])
{
  
  int i, n_elements, n_product;
  int *element;   // will contain a list of all elements in composition

  product_t     product;    // structure to hold product information
  composition_t propellant; // structure to hold propellant information

  double *sol;     // solution of the matrix
  double **matrix; // matrix to hold coefficient for the system resolution
  int    size;     // the size of the matrix


  int ch4, n2, hcn, h2, h2o;

  /* allocate memory to hold data */
  if (mem_alloc())
    return 1;
  
  // initialisation of the element list
  element = (int *)malloc(sizeof(int) * 100);
  for (i = 0; i < 100; i++)
    element[i] = -1;

  initialize_product(&product);

  load_thermo("thermo.dat");
  load_propellant("propellant.dat");
  

  /* Create a propellant */
  propellant.ncomp = 2;
  //propellant.molecule[0] = 34;  // AL
  //propellant.molecule[0] = 685; // O2
  //propellant.molecule[1] = 458; // H2
  propellant.molecule[0] = 766; // KCLO4
  propellant.molecule[1] = 788; // HTPB
  //propellant.coef[0] = GRAM_TO_MOL(10, propellant.molecule[0]);
  propellant.coef[0] = GRAM_TO_MOL(70, propellant.molecule[0]);
  propellant.coef[1] = GRAM_TO_MOL(30, propellant.molecule[1]);
  

  // n_elements contain the number of element
  n_elements = list_element(propellant, element); 
  printf("%d different elements in the propellant\n", n_elements);
  
  /* Print those elements */
  for (i = 0; i < n_elements; i++)
    printf("%s ", symb[element[i]] );
  printf("\n");
  

  // n_product contain the total number possible product of combustion
  // at temperature (1200)
  n_product = list_product(n_elements, element, &product, 1200);  

  printf("%d possible combustion product\n", n_product);
  printf("%d gazeous species\n", product.n[GAS]);
  printf("%d condensed species\n", product.n[CONDENSED]);
  

  // the size of the matrix to solve the problem
  size = n_elements + product.n[CONDENSED] + 1;

  // allocate the memory for the matrix
  matrix = (double **) calloc (size, sizeof(double *));
  for (i = 0; i < size; i++)
    matrix[i] = (double *) calloc (size+1, sizeof(double));

  sol = (double *) calloc (size, sizeof(double));


  fill_matrix(matrix, n_elements, element, &product, &propellant);
  print_matrix(matrix, size);
  gauss(matrix, sol, size);
  print_vec(sol, size);
 

  ch4 = thermo_search("C3H6");
  h2o = thermo_search("H2O(L)");
  n2  = thermo_search("N2");
  hcn = thermo_search("HCN");
  h2  = thermo_search("H");

  //print_thermo_list();
  //print_propellant_list();
    
  //print_thermo_info( h2 );


  //for (i = 0; i < 32; i++)
  //  printf("enthalpy à %d: %f\n", i*100, enthalpy( h2, i*100));


  //printf("enthalpy à 999: %f\n", enthalpy( ch4, 999));
  //printf("enthalpy à 1500: %f\n", enthalpy( ch4, 1500));
  //printf("enthalpy à 1500: %f\n", enthalpy( hcn, 1500));
  //printf("entropy: %f\n", entropy( ch4, 1500));
  //printf("specific heat: %f\n", specific_heat( ch4 , 300));


  //printf("%f\n",  potential(h2, 4000));



  /*
    temp = enthalpy(co2, 1700);
    printf("%f\n", delta_enthalpy(co2, 1800));
    printf("%f\n", enthalpy(co2, 1800));    
    printf("%f\n", entropy(co2, 1800));     
    printf("%f\n", specific_heat(co2, 1800));    
    printf("%f\n", gibbs(co2, 1800));


    printf("%f\n", heat_of_formation(840));
    printf("\n\nAverage time elapsed: %f\n", ((double)(finish - start)) 
    / (CLOCKS_PER_SEC * 20));
  */

  
  free(element);
  dealloc_product (&product);
  free (propellant_list);
  free (thermo_list);

  return 0;
  
}



