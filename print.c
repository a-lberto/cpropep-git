
#include <stdio.h>

#include "print.h"
#include "equilibrium.h"


extern propellant_t	*propellant_list;
extern thermo_t	    *thermo_list;

//extern float molar_mass[N_SYMB];

extern char symb[N_SYMB][3];


int print_propellant_info(int sp)
{
  int j;

  if (sp > MAX_PROPELLANT || sp < 0)
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

  if (sp > MAX_THERMO || sp < 0)
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
  for (i = 0; i < MAX_THERMO; i++)
    printf("%-4d %-15s % .2f\n", i, (thermo_list + i)->name,
           (thermo_list + i)->heat);
  
  return 0;
}

int print_propellant_list(void)
{
  int i;
  for (i = 0; i < MAX_PROPELLANT; i++)
    printf("%-4d %-30s %5d\n", i, (propellant_list + i)->name,
           (propellant_list +i)->heat);
 
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
  
  printf("%.4e mol of gaz\n", e->n);
  printf("molar fraction \t mol \t\t free energy\n");
  for (i = 0; i < e->p->n[GAS]; i++)
  {
    if (!(e->p->coef[GAS][i] == 0.0))
      printf("% .4e \t% .4e \t %f \t %s\n", 
             e->p->coef[GAS][i]/e->n,
             e->p->coef[GAS][i],
             gibbs_0(e->p->species[GAS][i], e->T),
             (thermo_list + e->p->species[GAS][i])->name);
  }
  if (e->p->n[CONDENSED] > 0)
    printf("Condensed species (mol)\n");
  for (i = 0; i < e->p->n[CONDENSED]; i++)
  {
    printf("%s  % .4e\n", (thermo_list + e->p->species[CONDENSED][i])->name,
           e->p->coef[CONDENSED][i]);
  }
  printf("Molar mass of product: % f g/mol\n", product_molar_mass(e));
  printf("Products enthalpy:     % f\n", product_enthalpy(e)*R*e->T);
  printf("Products entropy:      % f\n", product_entropy(e)*R);
         
  return 0;
}

int print_propellant_composition(equilibrium_t *e)
{
  int i, j;
  double enth = 0.0;
  
  printf("Propellant composition\n");
  printf("Code %-35s mol    Mass (g)  Composition\n", "Name");
  for (i = 0; i < e->c->ncomp; i++)
  {
    printf("%-4d  %-35s %.4f %.4f ", e->c->molecule[i],
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

  printf("Total mass: %f\n", propellant_mass(e));
  
  printf("\n");
  printf("Propellant properties\n");


  enth = propellant_enthalpy (e);

  /* not sure of the t, 273.15 or e->T */
  printf("Enthalpy: %.2f Joules  %.2f Joules/(RT)\n", enth*R*e->T, enth);
  
  printf("\n");
  return 0;
  
}


