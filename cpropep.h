#ifndef cpropep_h
#define cpropep_h

/* cpropep.h  -  Calculation of Complex Chemical Equilibrium           */
/* Copyright (C) 2000                                                  */
/* Antoine Lefebvre <antoine.lefebvre@polymtl.ca                       */
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


/**************************************************************
MACRO: Hold the number of species in for each type data type
***************************************************************/
#define MAX_PROPELLANT 1024
#define MAX_THERMO     1729

/* MACRO: Number of different molecule in propellant */
#define MAX_COMP       20

/* MACRO: Number of symbol in the symbol table */
#define N_SYMB         102

#define GRAM_TO_MOL(g, sp)   g/propellant_molar_mass(sp)

/****************************************************************
TYPE:  Enumeration of the possible state of a substance
*****************************************************************/
typedef enum 
{
  GAS,
  CONDENSED,
  STATE_LAST
} state_t;


/***************************************************************
TYPE: Structure to hold information of species contain in the
      thermo data file
****************************************************************/
typedef struct _thermo
{
  char    name[19];
  char    comments[57];
  int     nint;         /* number of different temperature interval */
  char    id[7];        /* identification code */
  int     elem[5]; 
  int     coef[5];
  state_t state;
  double  weight;       /* molecular weight */
  float   heat;         /* heat of formation at 298.15 K  (J/mol)  */
  double  dho;          /* HO(298.15) - HO(0) */
  float   range[4][2];  /* temperature range */
  int     ncoef[4];     /* number of coefficient for Cp0/R   */
  int     ex[4][8];     /* exponent in empirical equation */
  
  long double param[4][9];
  
  /* for species with data at only one temperature */
  /* especially condensed                          */
  float temp;
  float enth;
  
} thermo_t;

/***************************************************************
TYPE: Structure to hold information of species contain in the
      propellant data file
****************************************************************/
typedef struct _propellant
{
  char  name[30];  /* name of the propellant */
  int   elem[6];   /* element in the molecule (atomic number) max 6 */
  int   coef[6];   /* stochiometric coefficient of this element 
		      (0 for none) */
  int   heat; /* heat of formation in cal/gram */
  float density; /* density in pound/cubic inch */
  
} propellant_t;



/***************************************************************
TYPE: Hold the composition of a specific propellant
      ncomp is the number of component
      molecule[ ] hold the number in propellant_list corresponding
                  to the molecule
      coef[ ] hold the stochiometric coefficient

NOTE: It should be great to allocate the memory of the array in 
      function of the number of element

DATE: February 6, 2000
****************************************************************/
typedef struct _composition
{
  int ncomp;
  int molecule[MAX_COMP];
  float coef[MAX_COMP];
} composition_t;


/*****************************************************************
TYPE: Hold the composition of the combustion product. The molecule
      are separate between their different possible state.

NOTE: This structure should be initialize with the function 
      initialize_product.

DATE: February 13, 2000
******************************************************************/
typedef struct _product
{
  int    isalloc;              // true if the memory was allocated
  int    n[STATE_LAST];        // number of species for each state
  int    *species[STATE_LAST]; // list of species for each state
  double *coef[STATE_LAST];    // stoechiometric coefficient of each molecule
} product_t;



/***************************************************************
FUNCTION PROTOTYPE SECTION
****************************************************************/


/***************************************************************
FUNCTION: Print the information of a specie in the thermo_list

PARAMETER: an integer corresponding to the molecule

AUTHOR: Antoine Lefebvre
***************************************************************/
int print_thermo_info(int sp);

/*************************************************************
FUNCTION: Print the content of the respective list with the
          number which refer to the molecule

AUTHOR: Antoine Lefebvre
        modification bye Mark Pinese
**************************************************************/
int print_thermo_list(void);
int print_propellant_list(void);


/*************************************************************
FUNCTION: Search in the field name of thermo_list and return
          the value of the found item.

PARAMETER: A string corresponding to what we search, 
           example: "CO2"

COMMENTS: If nothing is found, it return -1

AUTHOR: Antoine Lefebvre
        modification bye Mark Pinese
**************************************************************/
int thermo_search(const char *str);



/************************************************************
FUNCTION: This function search for all elements present in
          the composition and fill the list with the 
	  corresponding number.

PARAMETER: comp is of type composition_t and hold the information
           about the propellant composition. list is a pointer
	   to an integer array.

COMMENTS: It is important that list is big enough to hold
          all the element, there is no check about how many
	  memory was allocated.

DATE: February 6, 2000

AUTHOR: Antoine Lefebvre
**************************************************************/
int list_element(composition_t comp, int *list);


/************************************************************
FUNCTION: This function search in thermo_list for all molecule
          that could be form with one or more of the element
	  in element_list. The function fill product_list with
	  the corresponding number of these molecule.

PARAMETER: n_element is the number of element in element_list,
           element_list is the list that was fill with list_element,
	   product_list is the list that will contain the result.
	   BE SURE PRODUCT_LIST IS BIG ENOUGH TO HOLD THE DATA,
	   there will be probably more than 200-300 molecule.

COMMENTS: The return value is the number of elements found

DATE: February 6, 2000

AUTHOR: Antoine Lefebvre
**************************************************************/
int list_product(int n_element, int *element_list, product_t *p, float T);


/*************************************************************
FUNCTION: Return the enthalpy of the molecule in thermo_list[sp]
          at the temperature T in K.

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It use the parametric form explain in the documentation
          to compute the value from the data read in the file
	  thermo.dat

AUTHOR: Antoine Lefebvre
**************************************************************/
double enthalpy(int sp, float T);

/*************************************************************
FUNCTION: Return the entropy of the molecule in thermo_list[sp]
          at the temperature T in K.

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It use the parametric form explain in the documentation
          to compute the value from the data read in the file
	  thermo.dat

AUTHOR: Antoine Lefebvre
**************************************************************/
double entropy(int sp, float T);

/*************************************************************
FUNCTION: Return the specific heat (Cp) of the molecule in 
          thermo_list[sp] at the temperature T in K.

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It use the parametric form explain in the documentation
          to compute the value from the data read in the file
	  thermo.dat

AUTHOR: Antoine Lefebvre
**************************************************************/
double specific_heat(int sp, float T);



/*************************************************************
FUNCTION: Return the variation of enthalpy of the molecule in 
          thermo_list[sp] between the temperature T in K and
	  298.15 K.

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It call enthalpy(...) for the enthalpy at temperature 
          T and use the field heat of thermo_t for the enthalpy
	  at 298.15

AUTHOR: Antoine Lefebvre
**************************************************************/
double delta_enthalpy(int sp, float T);

/*************************************************************
FUNCTION: Return the gibbs free energy of the molecule in 
          thermo_list[sp] at temperature T

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: g = H - ST where H is the enthalpy, T the temperature
          and S the entropy, so it call enthalpy(...) end
	  entropy(...)
**************************************************************/
double gibbs(int sb, state_t st, double nj, double n, float T, float P);



/* give the heat of formation of a propellant in kJ/mol */
double heat_of_formation(int molecule);

/* give the molar mass of a propellant or pruduct in g/mol */
double propellant_molar_mass(int molecule);
//double product_molar_mass(int molecule);


/**************************************************************
FUNCTION: This function initialize the the product structure.
          It allocate memory for the pointer and initialize
	  some value to -1

PARAMETER: A pointer to a structure product_t

AUTHOR: Antoine Lefebvre

DATE: February 13, 2000
****************************************************************/
int initialize_product(product_t *p);


/***************************************************************
FUNCTION: This function free all the pointer allocated in the
          product_t structure bye the initialisaztion
***************************************************************/
int dealloc_product(product_t *p);

#endif









