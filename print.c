
#include <stdio.h>

#include "print.h"
#include "performance.h"
#include "equilibrium.h"


extern propellant_t	*propellant_list;
extern thermo_t	    *thermo_list;

extern double g;

extern char symb[N_SYMB][3];

extern FILE * errorfile;
extern FILE * outputfile;

/* the minimum concentration we are interest to see */
#define CONC_MIN 1.0e-4 

#define ATM_TO_PSI 14.695949

int print_propellant_info(int sp)
{
  int j;

  if (sp > num_propellant || sp < 0)
    return -1;
  
  fprintf(outputfile, "Code %-35s Enthalpy  Density  Composition\n", "Name");
  fprintf(outputfile, "%d  %-35s % .2f % .2f", sp,
          (propellant_list + sp)->name,
          (float)(propellant_list + sp)->heat,
          (propellant_list + sp)->density);
  
  fprintf(outputfile, "  ");
  /* print the composition */
  for (j = 0; j < 6; j++)
  {
    if (!((propellant_list + sp)->coef[j] == 0))
      fprintf(outputfile, "%d%s ", (propellant_list + sp)->coef[j],
             symb[ (propellant_list + sp)->elem[j] ]);
  }
  fprintf(outputfile, "\n");
  return 0;
}

int print_thermo_info(int sp)
{
  int   i, j;
  thermo_t *s;

  if (sp > num_thermo || sp < 0)
    return -1;

  s = (thermo_list + sp);
  
  fprintf(outputfile, "---------------------------------------------\n");
  fprintf(outputfile, "Name: \t\t\t%s\n", s->name);
  fprintf(outputfile, "Comments: \t\t%s\n", s->comments);
  fprintf(outputfile, "Id: \t\t\t%s\n", s->id);
  fprintf(outputfile, "Chemical formula:\t");
  
  for (i = 0; i < 5; i++)
  {
    if (!(s->coef[i] == 0))
      fprintf(outputfile, "%d%s", s->coef[i], symb[ s->elem[i]]);
  }
  fprintf(outputfile, "\n");
  fprintf(outputfile, "State:\t\t\t");
  switch (s->state)
  {
    case GAS:
        fprintf(outputfile, "GAZ\n");
        break;
    case CONDENSED:
        fprintf(outputfile, "CONDENSED\n");
        break;
    default:
        printf("UNRECOGNIZE\n");
  }
  
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Molecular weight: \t\t% f g/mol\n", s->weight);
  fprintf(outputfile, "Heat of formation at 298.15 K : % f J/mol\n", s->heat);
  fprintf(outputfile, "HO(298.15) - HO(0): \t\t% f J/mol\n", s->dho);
  fprintf(outputfile, "Number of temperature range: % d\n\n", s->nint);
  
  for (i = 0; i < s->nint; i++)
  {
    fprintf(outputfile, "Interval: %f - %f \n", s->range[i][0],
            s->range[i][1]);
    for (j = 0; j < 9; j++)
      fprintf(outputfile, "% Le ", s->param[i][j]);
    fprintf(outputfile, "\n\n");
  }
  fprintf(outputfile, "---------------------------------------------\n");
  return 0;
}


int print_thermo_list(void)
{
  int i;
  for (i = 0; i < num_thermo; i++)
    fprintf(outputfile, "%-4d %-15s % .2f\n", i, (thermo_list + i)->name,
            (thermo_list + i)->heat);
  
  return 0;
}

int print_propellant_list(void)
{
  int i;
  for (i = 0; i < num_propellant; i++)
    fprintf(outputfile, "%-4d %-30s %5d\n", i, (propellant_list + i)->name,
            (propellant_list +i)->heat);
 
  return 0;
}


int print_condensed(product_t p)
{
  int i;
  for (i = 0; i < p.n[CONDENSED]; i ++)
    fprintf(outputfile, "%s ",
            (thermo_list + p.species[i][CONDENSED])->name );
  fprintf(outputfile, "\n");
  return 0;
}

int print_gazeous(product_t p)
{
  int i;
  for (i = 0; i < p.n[GAS]; i++)
    fprintf(outputfile, "%s ", (thermo_list + p.species[i][GAS])->name );
  fprintf(outputfile, "\n");
  return 0;
}

int print_product_composition(equilibrium_t *e)
{
  int i;
  
  fprintf(outputfile, "molar fraction \tmol \n");
  for (i = 0; i < e->p.n[GAS]; i++)
  {
    if (e->p.coef[i][GAS]/e->n > CONC_MIN)
    {
      fprintf(outputfile, "%.4e \t%.4e \t %s\n", 
              e->p.coef[i][GAS]/e->n,
              e->p.coef[i][GAS],
              (thermo_list + e->p.species[i][GAS])->name);
    }
  }
  
  if (e->p.n[CONDENSED] > 0)
    fprintf(outputfile, "Condensed species (mol)\n");
  for (i = 0; i < e->p.n[CONDENSED]; i++)
  {
    fprintf(outputfile, "         \t%.4e \t %s\n",
            e->p.coef[i][CONDENSED],
            (thermo_list + e->p.species[i][CONDENSED])->name);
  }
  fprintf(outputfile, "\n");
  return 0;
}

int print_product_properties(equilibrium_t *e)
{
  if (e->short_output)
  {
    fprintf(outputfile,
"T(K)   T(F)   P(atm) P(psi)  Enthalpy   Entropy Gas   RT/V   Molar mass\n");
    fprintf(outputfile,
            "%4.1f %4.1f %-6.2f %-7.2f %-10.2f %-7.2f %4.3f %-6.3f %4.3f\n",
            e->T, e->T*1.8-459.67, e->P, e->P*ATM_TO_PSI,
            product_enthalpy(e)*R*e->T,
            product_entropy(e)*R, e->n, e->P/e->n,
            product_molar_mass(e));
  }
  else
  {
    fprintf(outputfile, "Molar mass of product: % 11.3f g/mol\n",
            product_molar_mass(e));
    fprintf(outputfile, "Products enthalpy    : % 11.3f J\n",
            product_enthalpy(e)*R*e->T);
    fprintf(outputfile, "Products entropy     : % 11.3f J/K\n",
            product_entropy(e)*R);
    fprintf(outputfile, "Temperature          : % 11.3f K\n", e->T);
    fprintf(outputfile, "Pressure             : % 11.3f atm\n", e->P);
    fprintf(outputfile, "RT/V                 : % 11.3f\n", e->P/e->n);
    fprintf(outputfile, "Moles of gas         : % 11.3f\n", e->n);
  }
  
  return 0;
}

int print_propellant_composition(equilibrium_t *e)
{
  int i, j;
  
  fprintf(outputfile, "Propellant composition\n");
  fprintf(outputfile, "Code  %-35s mol    Mass (g)  Composition\n", "Name");
  for (i = 0; i < e->c.ncomp; i++)
  {
    fprintf(outputfile, "%-4d  %-35s %.4f %.4f ", e->c.molecule[i],
            (propellant_list + e->c.molecule[i])->name,
            e->c.coef[i], 
            e->c.coef[i]*propellant_molar_mass( e->c.molecule[i] ) );
    
    fprintf(outputfile, "  ");
    /* print the composition */
    for (j = 0; j < 6; j++)
    {
      if (!((propellant_list + e->c.molecule[i])->coef[j] == 0))
        fprintf(outputfile, "%d%s ",
                (propellant_list + e->c.molecule[i])->coef[j],
                symb[ (propellant_list + e->c.molecule[i])->elem[j] ]);
    }
    fprintf(outputfile, "\n");
  }

  fprintf(outputfile, "Total mass: % .2f g\n", propellant_mass(e));
  
  fprintf(outputfile, "Enthalpy  : % .2f Joules\n",
          propellant_enthalpy(e)*R*e->T);
  
  fprintf(outputfile, "\n");
  return 0;
  
}

int print_derivative_results(deriv_t *d)
{
  fprintf(outputfile, "Cp                   : % .2f J/Kg/K\n", d->cp);
  fprintf(outputfile, "Cv                   : % .2f J/Kg/K\n", d->cv);
  fprintf(outputfile, "Cp/Cv                : % .2f \n", d->cp_cv);
  fprintf(outputfile, "Isentropic exponent  : % .2f \n", d->isex);
  return 0;
}



int print_performance_information(performance_t *p)
{

  if (p->frozen_ok)
  {
    fprintf(outputfile,
            "\n--- Frozen equilibrium performance characteristics. ---\n");
    fprintf(outputfile,
         "         T(K)      P(atm)    Cp        Cp/Cv     Flow velocity\n");
    fprintf(outputfile,
            "Chamber: %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f\n",
            p->frozen.chamber.temperature,
            p->frozen.chamber.pressure,
            p->frozen.cp, p->frozen.cp_cv, 0.0);
    fprintf(outputfile,
            "Throat:  %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f\n",
            p->frozen.throat.temperature,
            p->frozen.throat.pressure,
            p->frozen.cp, p->frozen.cp_cv, p->frozen.throat.velocity);
    fprintf(outputfile,
            "Exit:    %-9.2f %-9.2f %-9.2f %-9.2f %-9.2f\n",
            p->frozen.exit.temperature,
            p->frozen.exit.pressure,
            p->frozen.cp, p->frozen.cp_cv, p->frozen.exit.velocity);
    
    fprintf(outputfile, "Specific impulse  : % 9.2f s\n",
            p->frozen.specific_impulse);

    fprintf(outputfile, "Exit to throat area ratio: %9.2f\n",
            (p->frozen.exit.temperature *
             p->frozen.throat.pressure *
             p->frozen.throat.velocity) /
            (p->frozen.throat.temperature *
             p->frozen.exit.pressure *
             p->frozen.exit.velocity ));
  }
  else
  {
    fprintf(outputfile, "Frozen performance was not computed\n");
  }

  if (p->equilibrium_ok)
  {
    fprintf(outputfile,
            "\n--- Shifting equilibrium performance characteristics. ---\n");
    fprintf(outputfile,
"         T(K)      P(atm) Cp        Cp/Cv  Is ex   Flow velocity  Molar mass\n");
    
    fprintf(outputfile,
            "Chamber: %-9.2f %-6.2f %-9.2f %-6.2f %-7.2f %-14.2f %-9.2f\n",
            p->equilibrium.chamber.temperature,
            p->equilibrium.chamber.pressure,
            p->equilibrium.chamber.deriv.cp,
            p->equilibrium.chamber.deriv.cp_cv,
            p->equilibrium.chamber.deriv.isex, 0.0,
            p->equilibrium.chamber.molar_mass);

    fprintf(outputfile,
            "Throat:  %-9.2f %-6.2f %-9.2f %-6.2f %-7.2f %-14.2f %-9.2f\n",
            p->equilibrium.throat.temperature,
            p->equilibrium.throat.pressure,
            p->equilibrium.throat.deriv.cp,
            p->equilibrium.throat.deriv.cp_cv,
            p->equilibrium.throat.deriv.isex,
            p->equilibrium.throat.velocity,
            p->equilibrium.throat.molar_mass);
            
    fprintf(outputfile,
            "Exit:    %-9.2f %-6.2f %-9.2f %-6.2f %-7.2f %-14.2f %-9.2f\n",
            p->equilibrium.exit.temperature,
            p->equilibrium.exit.pressure,
            p->equilibrium.exit.deriv.cp,
            p->equilibrium.exit.deriv.cp_cv,
            p->equilibrium.exit.deriv.isex,
            p->equilibrium.exit.velocity,
            p->equilibrium.exit.molar_mass);
    
    fprintf(outputfile, "Specific impulse  : % 9.2f s\n",
            p->equilibrium.specific_impulse);

    fprintf(outputfile, "Exit to throat area ratio: %9.2f\n\n",
            (p->equilibrium.exit.molar_mass *
             p->equilibrium.exit.temperature *
             p->equilibrium.throat.pressure *
             p->equilibrium.throat.velocity) /
            (p->equilibrium.throat.molar_mass *
             p->equilibrium.throat.temperature *
             p->equilibrium.exit.pressure *
             p->equilibrium.exit.velocity ));
  }
  else
  {
    fprintf(outputfile, "Equilibrium performance was not computed\n\n");
  }

  return 0;
}






