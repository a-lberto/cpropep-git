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

#include "libnum.h"


/* global variable containing the information about chemical species */
propellant_t	*propellant_list;
thermo_t	*thermo_list;


/***************************************************************************
Initial format of thermo.dat:
interval   variable   type	size	description
-----------------------------------------------------------------------------
(0, 18)    name	      string	18	compound name
(18, 73)   comments   string	55	comment
(73, 75)   nint	      int	2	the number of temperature intervals
(75, 81)   id	      string	6	the material id
81	   state      int	1	0 - GAS, else CONDENSED
(82, 95)   weight     float	13	molecular weight
(95, 108)  enth/heat  float	13	enthaply if nint == 0 
                                        else heat of formation
			...
			rest of file
			...
***************************************************************************/

int load_thermo(char *filename)
{
  FILE *fd;
  
  int i = 0;
  int j, k, l;
  
  char buf[88], *buf_ptr, tmp[32], *tmp_ptr;
  buf_ptr = &buf[0];
  tmp_ptr = &tmp[0];
  
  /* open the file for reading */
  if ((fd = fopen(filename, "r")) == NULL )
    return 1;
  
  
  while ((fgets(buf_ptr, 88, fd)) != NULL)
  {
    /* if the line is not commented */
    if (*(buf_ptr) != '!')
    {
      /* Read in the name and the comments */
      strncpy((thermo_list + i)->name, buf_ptr, 18);
      trim_spaces((thermo_list + i)->name, 18);
      
      strncpy((thermo_list + i)->comments, buf_ptr + 18, 55);
      trim_spaces((thermo_list + i)->comments, 55);
      
      if ((fgets(buf_ptr, 88, fd)) == NULL)  /* get a new line */
      {
        /* end of file occur */
	break;
      }
      
      strncpy(tmp_ptr, buf_ptr, 3);
      (thermo_list + i)->nint = atoi(tmp_ptr);
      
      strncpy((thermo_list + i)->id, buf_ptr + 3, 6);
      trim_spaces((thermo_list + i)->id, 6);
      
      /* get the chemical formula and coefficient */
      /* grep the elements ( 5  max )*/
      for (k = 0; k < 5; k++)
      {
	tmp[0] = buf[k * 8 + 10];
	tmp[1] = buf[k * 8 + 11];
	tmp[2] = '\0';
		    
	/* find the atomic number of the element */
	for (l = 0; l < N_SYMB; l++)
	{
	  if (!strcmp(tmp, symb[l]))
	  {
	    (thermo_list + i)->elem[k] = l;
	    break;
	  };
	};
		    
        // And the number of atoms
	strncpy(tmp_ptr, buf_ptr + k * 8 + 13, 6);
	tmp[6] = '\0';
	
	// Should this be an int?  If so, why is it stored in x.2 format?
	(thermo_list + i)->coef[k] = (int) atof(tmp_ptr);
      }
	       
      /* grep the state */
      if (buf[51] == '0')
        (thermo_list + i)->state = GAS;
      else
        (thermo_list + i)->state = CONDENSED;
      
      /* grep the molecular weight */
      strncpy(tmp_ptr, buf_ptr + 52, 13);
      tmp[13] = '\0';
      (thermo_list + i)->weight = atof(tmp_ptr);
      
      /* grep the heat of formation (J/mol) or enthalpy if condensed */
      /* The values are assigned in the if block following */
      strncpy(tmp_ptr, buf_ptr + 65, 13);
      tmp[13] = '\0';
      
      /* now get the data */
      /* there is '(thermo_list + i)->nint' set of data */
      if ((thermo_list + i)->nint == 0)
      {
	/* Set the enthalpy */
	(thermo_list + i)->enth = atof(tmp_ptr);
	
	/* condensed phase, different info */
	if ((fgets(buf_ptr, 88, fd)) == NULL) 
	{
	  /* end of file occur */
	  break;
	}
			  
	/* treat the line */
	/* get the temperature of the assigned enthalpy */
	strncpy(tmp_ptr, buf_ptr + 1, 10);
	tmp[10] = '\0';
	
	(thermo_list + i)->temp = atof(tmp_ptr);
      }
      else 
      { 
	/* Set the heat of formation */
	(thermo_list + i)->heat = atof(tmp_ptr);
	
	for (j = 0; j < (thermo_list + i)->nint; j++)
	{
	  /* Get the first line of three */
	  if ( (fgets(buf_ptr, 88, fd)) == NULL) 
	  {
	    /* end of file occur */
	    break;
	  }
	  
	  /* low */
	  strncpy(tmp_ptr, buf_ptr + 1, 10);
	  tmp[10] = '\0';
	  (thermo_list + i)->range[j][0] = atof(tmp_ptr);
	  
	  /* high */
	  strncpy(tmp_ptr, buf_ptr + 11, 10);
	  tmp[10] = '\0';
	  (thermo_list + i)->range[j][1] = atof(tmp_ptr);
	  
	  tmp[0] = buf[22];
	  tmp[1] = '\0';
	  (thermo_list + i)->ncoef[j] = atoi(tmp_ptr);
	  
	  /* grep the exponent */
	  for (l = 0; l < 8; l++)
	  {
	    strncpy(tmp_ptr, buf_ptr + l * 5 + 23, 5);
	    tmp[5] = '\0';					     
	    (thermo_list + i)->ex[j][l] = atoi(tmp_ptr);
	  }
	  
	  /* HO(298.15) -HO(0) */
	  strncpy(tmp_ptr, buf_ptr + 65, 15);
	  tmp[15] = '\0';
	  (thermo_list + i)->dho = atof(tmp);
	  
	  
	  /* Get the second line of three */
	  if ( (fgets(buf_ptr, 88, fd)) == NULL) 
	  {
	    /* end of file occur */
	    break;
	  }
			       
	  /* grep the first data line */
	  /* there are 5 coefficients */
	  for (l = 0; l < 5; l++)
	  {
	    strncpy(tmp_ptr, buf_ptr + l * 16, 16);
	    tmp[16] = '\0';
	    
	    (thermo_list + i)->param[j][l] = atof(tmp_ptr);
	  }
	  
	  /* Get the third line of three */
	  if ( (fgets(buf_ptr, 88, fd)) == NULL) 
	  {
	    /* end of file occur */
	    break;
	  }
	  
	  /* grep the second data line */
	  for (l = 0; l < 2; l++)
	  {
	    strncpy(tmp_ptr, buf_ptr + l * 16, 16);
	    tmp[16] = '\0';
	    
	    (thermo_list + i)->param[j][l + 5] = atof(tmp_ptr);
	    
	  }
	  
	  for (l = 0; l < 2; l++)
	  {
	    strncpy(tmp_ptr, buf_ptr + l * 16 + 48, 16);
	    tmp[16] = '\0';
	    
	    (thermo_list + i)->param[j][l + 7] = atof(tmp_ptr);
	    
	  }
	}
      }
      i++;
    }
  }
  
  fclose(fd);
  return i;

}


int load_propellant(char *filename) 
{
  
  FILE *fd;
  
  int i = 0;
  int j, k;
  
  /* temporary string to store string in order to treat the informations */
  char buf[88], *buf_ptr, tmp[32], *tmp_ptr;
  buf_ptr = &buf[0];
  tmp_ptr = &tmp[0];
  
  /* open the file for reading */
  if ((fd = fopen(filename, "r")) == NULL )
    return 1;
  
  while ((fgets(buf_ptr, 88, fd)) != NULL)
  {
    /* if the line is not commented */
    if (!((*(buf_ptr) == '*') || (*(buf_ptr) == '+')))
    {  
      /* grep the name */
      strncpy((propellant_list + i)->name, buf_ptr + 9, 29);
      trim_spaces((propellant_list + i)->name, 29);
      
      for (j = 0; j < 6; j++)
      {
	tmp[0] = buf[j * 5 + 39];
	tmp[1] = buf[j * 5 + 40];
	tmp[2] = buf[j * 5 + 41];
	tmp[3] = '\0';
	
	(propellant_list + i)->coef[j] = atoi(tmp);
	
	tmp[0] = buf[j * 5 + 42];
	tmp[1] = buf[j * 5 + 43];
	tmp[2] = '\0';
		   
	      /* find the atomic number of the element */
	for (k = 0; k < N_SYMB; k++)
	{
	  if (!(strcmp(tmp, symb[k]))) 
	  {
	    (propellant_list + i)->elem[j] = k;
	    break;
	  }
	}
      }
	  
      strncpy(tmp_ptr, buf_ptr + 70, 5);
      tmp[5] = '\0';		    
      propellant_list[i].heat = atoi(tmp);
      
      strncpy(tmp_ptr, buf_ptr + 70, 5);
      tmp[5] = '\0';
      propellant_list[i].density = atof(tmp);
      
      i++;
    }
  }  
  fclose(fd);     
  return i;
}


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
      val = -s.param[i][0] * pow(T, -2) + s.param[i][1] * pow(T, -1) * log(T)
	+ s.param[i][2] + s.param[i][3] * T / 2 + s.param[i][4] * pow(T, 2) / 3
	+ s.param[i][5] * pow(T, 3) / 4 + s.param[i][6] * pow(T, 4) / 5 + s.param[i][7] / T;
      
      return val * R * T; /* to convert from dimensionless to J/mol */
    }
  }
  
  printf("Error: temperature out of range\n");
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
  
  printf("Error: temperature out of range\n");
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
  
  printf("Error: temperature out of range\n");
  return 0;
}

double delta_enthalpy(int sp, float T)
{
  /* delta henthalpy in J/mol */
  return enthalpy(sp, T) - (thermo_list + sp)->heat;
}

double gibbs(int sp, float T)
{
  /* G = H -TS */ 
  return enthalpy(sp, T) - T * entropy(sp, T);
}

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


void trim_spaces(char *str, unsigned int len)
{
  /* Removes trailing ' ' characters.  If the entire string is ' ' characters, leaves one. */
  unsigned int i;
  
  for (i = len - 1; i > 0; i--)
  {
    if (*(str + i) != ' ')
    {
      *(str + i + 1) = '\0';
      return;
    }
  }
  
  *(str + 1) = '\0';
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

int list_product(int n_element, int *element_list, product_t *p)
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
      p->species[st][ p->n[st] ] = j;
      p->n[st]++;
      e++;
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

  // initialize the coef to the first approximation
  for (j = 0; j < STATE_LAST; j++)
  {
    for (i = 0; i < p->n[j]; i++)
      p->coef[j][i] = 0.1/p->n[j];
  }
  
  return e;
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

int element_coef(int element, int molecule)
{
  int i;
  for (i = 0; i < 5; i++)
  {
    if (thermo_list[molecule].elem[i] == element)
      return thermo_list[molecule].coef[i];
  }
  return 0;
}


int main(int argc, char *argv[])
{
  /* clock_t start, finish; */
  
  int i, j, k, n_elements, n_product;

  double tmp;

  int *element;   // will contain a list of all elements in composition

  product_t product;        // structure to hold product information
  composition_t propellant; // structure to hold propellant information

  double **matrix; // matrix to hold coefficient for the system resolution
  int size;


  //double temp;
  //int co2;
  
  /* allocate memory to hold data */
  if (mem_alloc())
    return 1;
  
  // initialisation of the element list
  element = (int *)calloc(100, sizeof(int));
  for (i = 0; i < 100; i++)
    element[i] = -1;

  initialize_product(&product);

 
      
  load_thermo("thermo.dat");
  load_propellant("propellant.dat");
  


  propellant.ncomp = 2;
  
  //propellant.molecule[0] = 34;  // AL
  
  propellant.molecule[0] = 685; // O2
  propellant.molecule[1] = 458; // H2

  //propellant.molecule[0] = 766; // KCLO4
  //propellant.molecule[1] = 788; // HTPB
  
  //propellant.coef[0] = GRAM_TO_MOL(10, propellant.molecule[0]);
  propellant.coef[0] = GRAM_TO_MOL(70, propellant.molecule[0]);
  propellant.coef[1] = GRAM_TO_MOL(30, propellant.molecule[1]);
  
  
  n_elements = list_element(propellant, element); // ele contain the number of element
  printf("%d\n", n_elements);
  
  for (i = 0; i < n_elements; i++)
    printf("%s ", symb[element[i]] );
  
  printf("\n");
  
  n_product = list_product(n_elements, element, &product);  // pro contain the total number
  // of product
  printf("%d\n", n_product);
  printf("%d\n", product.n[GAS]);
  printf("%d\n", product.n[CONDENSED]);
  
  size = n_elements + product.n[CONDENSED] + 2;

  // allocate the memory for the matrix
  matrix = (double **) malloc (sizeof(double *) * size);
  for (i = 0; i < size; i++)
    matrix[i] = (double *) malloc (sizeof(double) * (size+1));

  for (i = 0; i < n_elements; i++) // each column
  {
    for (j = 0; j < n_elements; j++) // each row
    {
      tmp = 0.0;
      for (k = 0; k < product.n[GAS]; k++)
  	tmp += element_coef( element[j], product.species[GAS][k]) * 
	  element_coef( element[i], product.species[GAS][k] ) * 
	  product.coef[GAS][k];
      
      matrix[j][i] = tmp;

    }
  }

  print_square_matrix(matrix, n_elements);
 



  
  //co2 = thermo_search("CO2");
    
  //print_thermo_list();
  //  print_propellant_list();
    /*
    print_thermo_info( co2 );
    
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



