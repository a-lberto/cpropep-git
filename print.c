
#include <stdio.h>

#include "print.h"
#include "equilibrium.h"


extern propellant_t	*propellant_list;
extern thermo_t	    *thermo_list;

//extern float molar_mass[N_SYMB];

extern char symb[N_SYMB][3];

/* the minimum concentration we are interest to see */
#define CONC_MIN 1e-4 

int print_propellant_info(int sp)
{
  int j;

  if (sp > num_propellant || sp < 0)
    return -1;
  
  printf("Code %-35s Enthalpy  Density  Composition\n", "Name");
  printf("%d  %-35s % .2f % .2f", sp,
         (propellant_list + sp)->name,
         (float)(propellant_list + sp)->heat,
         (propellant_list + sp)->density);
  
  printf("  ");
  /* print the composition */
  for (j = 0; j < 6; j++)
  {
    if (!((propellant_list + sp)->coef[j] == 0))
      printf("%d%s ", (propellant_list + sp)->coef[j],
             symb[ (propellant_list + sp)->elem[j] ]);
  }
  printf("\n");
  return 0;
}

int print_thermo_info(int sp)
{
  int   i, j;
  thermo_t *s;

  if (sp > num_thermo || sp < 0)
    return -1;

  s = (thermo_list + sp);
  
  printf("---------------------------------------------\n");
  printf("Name: \t\t\t%s\n", s->name);
  printf("Comments: \t\t%s\n", s->comments);
  printf("Id: \t\t\t%s\n", s->id);
  printf("Chemical formula:\t");
  
  for (i = 0; i < 5; i++)
  {
    if (!(s->coef[i] == 0))
      printf("%d%s", s->coef[i], symb[ s->elem[i]]);
  }
  printf("\n");
  printf("State:\t\t\t");
  switch (s->state)
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
  printf("Molecular weight: \t\t% f g/mol\n", s->weight);
  printf("Heat of formation at 298.15 K : % f J/mol\n", s->heat);
  printf("HO(298.15) - HO(0): \t\t% f J/mol\n", s->dho);
  printf("Number of temperature range: % d\n\n", s->nint);
  
  for (i = 0; i < s->nint; i++)
  {
    printf("Interval: %f - %f \n", s->range[i][0], s->range[i][1]);
    for (j = 0; j < 9; j++)
      printf("% Le ", s->param[i][j]);
    printf("\n\n");
  }
  printf("---------------------------------------------\n");
  return 0;
}


int print_thermo_list(void)
{
  int i;
  for (i = 0; i < num_thermo; i++)
    printf("%-4d %-15s % .2f\n", i, (thermo_list + i)->name,
           (thermo_list + i)->heat);
  
  return 0;
}

int print_propellant_list(void)
{
  int i;
  for (i = 0; i < num_propellant; i++)
    printf("%-4d %-30s %5d\n", i, (propellant_list + i)->name,
           (propellant_list +i)->heat);
 
  return 0;
}


int print_condensed(product_t p)
{
  int i;
  for (i = 0; i < p.n[CONDENSED]; i ++)
    printf("%s ",  (thermo_list + p.species[i][CONDENSED])->name );
  printf("\n");
  return 0;
}

int print_gazeous(product_t p)
{
  int i;
  for (i = 0; i < p.n[GAS]; i++)
    printf("%s ", (thermo_list + p.species[i][GAS])->name );
  printf("\n");
  return 0;
}

int print_product_composition(equilibrium_t *e)
{
  int i;
  
  printf("%.4e mol of gaz\n", e->n);
  printf("molar fraction \t mol \n");
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    if (e->p.coef[GAS][i]/e->n > CONC_MIN)
      printf("% .4e \t% .4e \t %s\n", 
             e->p.coef[i][GAS]/e->n,
             e->p.coef[i][GAS],
             (thermo_list + e->p.species[i][GAS])->name);
  }
  if (e->p.n[CONDENSED] > 0)
    printf("Condensed species (mol)\n");
  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {
    printf("%s  % .4e\n", (thermo_list + e->p.species[i][CONDENSED])->name,
           e->p.coef[i][CONDENSED]);
  }
  printf("\n");
  return 0;
}

int print_product_properties(equilibrium_t *e)
{
  printf("Molar mass of product: % .2f g/mol\n", product_molar_mass(e));
  printf("Products enthalpy    : % .2f J\n",     product_enthalpy(e)*R*e->T);
  printf("Products entropy     : % .2f J/K\n",   product_entropy(e)*R);
  printf("Temperature          : % .2f K\n",     e->T);
  printf("Pressure             : % .2f atm\n",   e->P);
  return 0;
}

int print_propellant_composition(equilibrium_t *e)
{
  int i, j;
  double enth = 0.0;
  
  printf("Propellant composition\n");
  printf("Code %-35s mol    Mass (g)  Composition\n", "Name");
  for (i = 0; i < e->c.ncomp; i++)
  {
    printf("%-4d  %-35s %.4f %.4f ", e->c.molecule[i],
           (propellant_list + e->c.molecule[i])->name,
           e->c.coef[i], 
           e->c.coef[i]*propellant_molar_mass( e->c.molecule[i] ) );
    
    printf("  ");
    /* print the composition */
    for (j = 0; j < 6; j++)
    {
      if (!((propellant_list + e->c.molecule[i])->coef[j] == 0))
        printf("%d%s ", (propellant_list + e->c.molecule[i])->coef[j],
               symb[ (propellant_list + e->c.molecule[i])->elem[j] ]);
    }
    printf("\n");
  }

  printf("Total mass: % .2f g\n", propellant_mass(e));
  
  printf("\n");
  printf("Propellant properties\n");
  enth = propellant_enthalpy (e);
  printf("Enthalpy: %.2f Joules\n", enth*R*e->T);
  
  printf("\n");
  return 0;
  
}

int print_derivative_results(deriv_t d)
{
  printf("Cp                   : % .2f \n", d.cp);
  printf("Cv                   : % .2f \n", d.cv);
  printf("Cp/Cv                : % .2f \n", d.cp_cv);
  printf("Isentropic exponent  : % .2f \n", d.isex);
  return 0;
}



