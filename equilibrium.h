#ifndef equilibrium_h
#define equilibrium_h

/* equilibrium.h  -  Calculation of Complex Chemical Equilibrium           */
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
#define MAX_PROPELLANT 1030
#define MAX_THERMO     1729

/* MACRO: Number of different molecule in propellant */
#define MAX_COMP       20

/* MACRO: Number of symbol in the symbol table */
#define N_SYMB         102

#define GRAM_TO_MOL(g, sp)   g/propellant_molar_mass(sp)

#define __min(a, b) ( (a) <= (b) ? (a) : (b))
#define __max(a, b) ( (a) >= (b) ? (a) : (b))

#define _min(a, b, c) __min( __min(a, b), c)
#define _max(a, b, c) __max( __max(a, b), c)

extern int global_verbose;
extern const float R;

/****************************************************************
TYPE:  Enumeration of the possible state of a substance
*****************************************************************/
typedef enum 
{
  GAS,
  CONDENSED,
  STATE_LAST
} state_t;


typedef enum
{
  TP,          // assign temperature and pressure
  HP,          // assign enthalpy and pressure
  SP           // assign entropy and pressure
} problem_t;


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
  int   heat;      /* heat of formation in Joule/gram */
  float density;   /* density in g/cubic cm */
  
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
  int    ncomp;
  int    molecule[MAX_COMP];
  double coef[MAX_COMP];
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
  int    *species[STATE_LAST]; // list of possible species for each state
  double *coef[STATE_LAST];    // stoechiometric coefficient of each molecule

  int    *condensed_ok;        // list containing the condensed species
                               // that respect the criterium to be present
} product_t;


/*****************************************************************
TYPE: Hold important information relative to an equilibrium 
      calculation

NOTE: 

DATE: February 24, 2000
******************************************************************/
typedef struct _equilibrium
{
  unsigned int  verbose;      /* verbose level */
    
  composition_t *c;           /* pointer to a propellant composition */
  product_t     *p;           /* pointer to a product struct */

  /* list of element in the composition */
  int  *element;
  int   n_element;

  int is_state_set; /* true if you have set the state */
  double T; /* temperature */
  double P; /* pressure */

    /* not needed */
  int    isequil;             /* true when the equilibrium have been compute */

  double         n;           /* total number of mole */
  double         delta_ln_n;
  double        *delta_ln_nj; /* hold delta ln(nj) for each gazeous species */
} equilibrium_t;




/***************************************************************
FUNCTION PROTOTYPE SECTION
****************************************************************/

/**************************************************************
FUNCTION: This function allocate the memory for the two global
          list thermo_list and propellant_list

COMMENTS: It should be call before loading the thermo and
          the propellant file

AUTHOR: Antoine Lefebvre
***************************************************************/
int mem_alloc(void);



int set_verbose(equilibrium_t *e, int v);

/*************************************************************
FUNCTION: Search in the field name of thermo_list and return
          the value of the found item.

PARAMETER: A string corresponding to what we search, 
           example: "CO2"

COMMENTS: If nothing is found, it return -1

AUTHOR: Antoine Lefebvre
        modification bye Mark Pinese
**************************************************************/
int thermo_search(char *str);

int propellant_search(char *str);

/************************************************************
FUNCTION: This function search for all elements present in
          the composition and fill the list with the 
	  corresponding number.

PARAMETER: e is of type equilibrium_t and hold the information
           about the propellant composition.

COMMENTS: It fill the member element in equilibrium_t

DATE: February 6, 2000

AUTHOR: Antoine Lefebvre
**************************************************************/
int list_element(equilibrium_t *e);


/************************************************************
FUNCTION: This function search in thermo_list for all molecule
          that could be form with one or more of the element
	  in element_list. The function fill product_list with
	  the corresponding number of these molecule.

PARAMETER: e is a pointer to an equilibrium_t structure

COMMENTS: The return value is the number of elements found

DATE: February 6, 2000

AUTHOR: Antoine Lefebvre
**************************************************************/
int list_product(equilibrium_t *e);


/*************************************************************
FUNCTION: Return the enthalpy of the molecule in thermo_list[sp]
          at the temperature T in K. (Ho/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It use the parametric form explain in the documentation
          to compute the value from the data read in the file
	  thermo.dat

AUTHOR: Antoine Lefebvre
**************************************************************/
double enthalpy_0(int sp, float T);

/*************************************************************
FUNCTION: Return the entropy of the molecule in thermo_list[sp]
          at the temperature T in K. (So/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It use the parametric form explain in the documentation
          to compute the value from the data read in the file
	  thermo.dat

AUTHOR: Antoine Lefebvre
**************************************************************/
double entropy_0(int sp, float T);

double entropy(int sp, state_t st, double nj, double n, float T, float P);

/*************************************************************
FUNCTION: Return the specific heat (Cp) of the molecule in 
          thermo_list[sp] at the temperature T in K. (Cp/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: It use the parametric form explain in the documentation
          to compute the value from the data read in the file
	  thermo.dat

AUTHOR: Antoine Lefebvre
**************************************************************/
double specific_heat_0(int sp, float T);


/*************************************************************
FUNCTION: Return true if the thermochemical data are define for
          this temperature.

PARAMETER: The same as for entropy

COMMENTS:  It is useful to determine if a specie is present at
           a given temperature.

AUTHOR: Antoine Lefebvre
**************************************************************/
int temperature_check(int sp, float T);


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

double propellant_enthalpy(equilibrium_t *e);
double product_enthalpy(equilibrium_t *e);
double product_entropy(equilibrium_t *e);

double propellant_mass(equilibrium_t *e);

/*************************************************************
FUNCTION: Return the gibbs free energy of the molecule in 
          thermo_list[sp] at temperature T. (uo/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: g = H - ST where H is the enthalpy, T the temperature
          and S the entropy.
**************************************************************/
double gibbs_0(int sp, float T);


/*************************************************************
FUNCTION: Return the gibbs free energy of the molecule in 
          thermo_list[sp] at temperature T, pressure P. (u/RT)

PARAMETER: sp is the position in the array of the molecule
           T is the temperature in K

COMMENTS: g = uo + ln(nj/n) + ln(P) for gazes
          g = uo for condensed

AUTHOR: Antoine Lefebvre
**************************************************************/
double gibbs(int sp, state_t st, double nj, double n, float T, float P);


/***************************************************************
FUNCTION: Return the heat of formation of a propellant in kJ/mol
****************************************************************/
double heat_of_formation(int molecule);

/*************************************************************
FUNCTION: Return the molar mass of a propellant (g/mol)

PARAMETER: molecule is the number in propellant_list
**************************************************************/
double propellant_molar_mass(int molecule);


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
FUNCTION: This function initialize the equilibrium structure.
          The function allocate memory for all the structure
	  it need. It is important to call dealloc_equilibrium
	  after.

AUTHOR:   Antoine Lefebvre

DATE: February 27, 2000
****************************************************************/
int initialize_equilibrium(equilibrium_t *e);


/***************************************************************
FUNCTION: Dealloc what have been allocated by 
          initialize_equilibrium
***************************************************************/
int dealloc_equillibrium(equilibrium_t *e);


/***************************************************************
FUNCTION: This function free all the pointer allocated in the
          product_t structure bye the initialisaztion
***************************************************************/
int dealloc_product(product_t *p);

/***************************************************************
FUNCTION: Set the state at which we want to compute the 
          equilibrium.

PARAMETER: e is a pointer to an equilibrium_t structure
           T is the temperature in deg K
	   P is the pressure in atm

AUTHOR:    Antoine Lefebvre
****************************************************************/
int set_state(equilibrium_t *e, double T, double P);


/***************************************************************
FUNCTION: Add a new molecule in the propellant

PARAMETER: e is a pointer to the equilibrium_t structure
           sp is the number of the molecule in the list
	   mol is the quantity in mol

AUTHOR:    Antoine Lefebvre
****************************************************************/
int add_in_propellant(equilibrium_t *e, int sp, double mol);

/***************************************************************
FUNCTION: Return the stochiometric coefficient of an element
          in a molecule. If the element isn't present, it return 0.

COMMENTS: There is a different function for the product and for the
          propellant.

AUTHOR:   Antoine Lefebvre
****************************************************************/
int product_element_coef(int element, int molecule);
int propellant_element_coef(int element, int molecule);




/***************************************************************
FUNCTION: This function fill the matrix in function of the data
          store in the structure equilibrium_t. The solution
	  of this matrix give corresction to initial estimate.

COMMENTS: It use the theory explain in 
          "Computer Program for Calculation of Complex Chemical
	  Equilibrium Compositions, Rocket Performance, Incident
	  and Reflected Shocks, and Chapman-Jouguet Detonations"
	  by Gordon and McBride

AUTHOR:   Antoine Lefebvre
****************************************************************/
int fill_matrix(double **matrix, equilibrium_t *e, problem_t P);


/****************************************************************
FUNCTION: This function compute the equilibrium composition at
          at specific pressure/temperature point. It use fill_matrix
	  to obtain correction to initial estimate. It correct the 
	  value until equilibrium is obtain.

AUTHOR:   Antoine Lefebvre
******************************************************************/
int equilibrium(equilibrium_t *equil, problem_t P);


double product_molar_mass(equilibrium_t *e);

#endif









