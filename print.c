
#include <stdio.h>

#include "print.h"
#include "performance.h"
#include "equilibrium.h"
#include "conversion.h"
#include "thermo.h"
#include "const.h"


FILE * errorfile;
FILE * outputfile;

int print_propellant_info(int sp)
{
  int j;

  if (sp > num_propellant || sp < 0)
    return -1;
  
  fprintf(outputfile, "Code %-35s Enthalpy  Density  Composition\n", "Name");
  fprintf(outputfile, "%d  %-35s % .4f % .2f", sp,
          (propellant_list + sp)->name,
          (propellant_list + sp)->heat,
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
  fprintf(outputfile, "Assign enthalpy               : % f J/mol\n", s->enth);
  fprintf(outputfile, "HO(298.15) - HO(0): \t\t% f J/mol\n", s->dho);
  fprintf(outputfile, "Number of temperature range: % d\n\n", s->nint);
  
  for (i = 0; i < s->nint; i++)
  {
    fprintf(outputfile, "Interval: %f - %f \n", s->range[i][0],
            s->range[i][1]);
    for (j = 0; j < 9; j++)
      fprintf(outputfile, "% .9e ", s->param[i][j]);
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
    fprintf(outputfile, "%-4d %-30s %5f\n", i, (propellant_list + i)->name,
            (propellant_list +i)->heat);
 
  return 0;
}


int print_condensed(product_t p)
{
  int i;
  for (i = 0; i < p.n[CONDENSED]; i ++)
    fprintf(outputfile, "%s ",
            (thermo_list + p.species[CONDENSED][i])->name );
  fprintf(outputfile, "\n");
  return 0;
}

int print_gazeous(product_t p)
{
  int i;
  for (i = 0; i < p.n[GAS]; i++)
    fprintf(outputfile, "%s ", (thermo_list + p.species[GAS][i])->name );
  fprintf(outputfile, "\n");
  return 0;
}

int print_product_composition(equilibrium_t *e, short npt)
{
  int i, j, k;
  
  double mol_g = e->itn.n;
  
  for (i = 0; i < e->product.n[CONDENSED]; i++)
    mol_g += e->product.coef[CONDENSED][i];
  
  fprintf(outputfile, "Molar fractions\n\n");
  for (i = 0; i < e->product.n[GAS]; i++)
  {
    if (e->product.coef[GAS][i]/e->itn.n > 0.0)
    {
      fprintf(outputfile, "%-20s",
              (thermo_list + e->product.species[GAS][i])->name);
//              e->product.coef[GAS][i]/mol_g);

      for (j = 0; j < npt; j++)
        fprintf(outputfile, " %11.4e", (e+j)->product.coef[GAS][i]/mol_g);
      fprintf(outputfile,"\n");
      
    }
  }

  if (e->product.n[CONDENSED] > 0)
  {
    fprintf(outputfile, "Condensed species\n");
    for (i = 0; i < e->product.n[CONDENSED]; i++)
    {
      fprintf(outputfile,   "%-20s %.4e\n",
              (thermo_list + e->product.species[CONDENSED][i])->name,
              e->product.coef[CONDENSED][i]/mol_g);
    }
  }
  fprintf(outputfile, "\n");
  return 0;
}


int print_propellant_composition(equilibrium_t *e)
{
  int i, j;
  
  fprintf(outputfile, "Propellant composition\n");
  fprintf(outputfile, "Code  %-35s mol    Mass (g)  Composition\n", "Name");
  for (i = 0; i < e->propellant.ncomp; i++)
  {
    fprintf(outputfile, "%-4d  %-35s %.4f %.4f ", e->propellant.molecule[i],
            (propellant_list + e->propellant.molecule[i])->name,
            e->propellant.coef[i], 
            e->propellant.coef[i]*propellant_molar_mass(e->propellant.molecule[i]));
    
    fprintf(outputfile, "  ");
    /* print the composition */
    for (j = 0; j < 6; j++)
    {
      if (!((propellant_list + e->propellant.molecule[i])->coef[j] == 0))
        fprintf(outputfile, "%d%s ",
                (propellant_list + e->propellant.molecule[i])->coef[j],
                symb[ (propellant_list + e->propellant.molecule[i])->elem[j] ]);
    }
    fprintf(outputfile, "\n");
  }

  if (e->product.element_listed)
  {
    fprintf(outputfile, "%d different elements\n", e->product.n_element);
    /* Print those elements */
    for (i = 0; i < e->product.n_element; i++)
      fprintf(outputfile, "%s ", symb[e->product.element[i]] );
    fprintf(outputfile, "\n");
  }
  
  fprintf(outputfile, "Total mass: % f g\n", propellant_mass(e));
  
  fprintf(outputfile, "Enthalpy  : % .2f Joules\n",
          propellant_enthalpy(e)*propellant_mass(e));
  
  fprintf(outputfile, "\n");

  if (e->product.product_listed)
  {
    fprintf(outputfile, "%d possible gazeous species\n", e->product.n[GAS]);
    if (global_verbose > 1)
      print_gazeous(e->product);
    fprintf(outputfile, "%d possible condensed species\n", e->product.n[CONDENSED]);
    if (global_verbose > 1)
      print_condensed(e->product);
  }
  
  return 0;
}

int print_performance_information(equilibrium_t *e, short npt)
{
  short i;
  
  fprintf(outputfile, "Ae/At            :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.ae_at);
  fprintf(outputfile, "\n");
  
  fprintf(outputfile, "A/dotm           :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.a_dotm);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "C* (m/s)         :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.cstar);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "Cf               :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.cf);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "Ivac (m/s)       :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.Ivac);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "Isp (m/s)        :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.Isp);
  fprintf(outputfile, "\n");

  fprintf(outputfile, "Isg/g (s)        :            ");
  for (i = 1; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->performance.Isp/Ge);
  fprintf(outputfile, "\n");

  return 0;
}


int print_product_properties(equilibrium_t *e, short npt)
{
  short i;
  
  fprintf(outputfile, "Pressure (atm)   :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.P);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Temperature (K)  :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.T);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "H (kJ/kg)        :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.H);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "U (kJ/kg)        :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.U);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "G (kJ/kg)        :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.G);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "S (kJ/(kg)(K)    :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.S);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "M (g/mol)        :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.3f", (e+i)->properties.M);
  fprintf(outputfile, "\n");
  
  fprintf(outputfile, "(dLV/dLP)t       :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.dV_P);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "(dLV/dLT)p       :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.dV_T);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Cp (kJ/(kg)(K))  :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Cp);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Cv (kJ/(kg)(K))  :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Cv);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Cp/Cv            :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Cp/(e+i)->properties.Cv);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Gamma            :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Isex);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "Vson (m/s)       :");
  for (i = 0; i < npt; i++)
    fprintf(outputfile, " % 11.5f", (e+i)->properties.Vson);
  fprintf(outputfile, "\n");
  fprintf(outputfile, "\n");
  return 0;
}

/*
int print_performance_information(performance_t *p)
{
  
  if (p->frozen_ok)
  {
    fprintf(outputfile,
            "\n--- Frozen equilibrium performance characteristics. ---\n");

    
    printf("Ae/At      : %8.3f \t %8.3f\n", 1.0,
           (p->frozen.exit.temperature *
            p->frozen.throat.pressure *
            p->frozen.throat.velocity) /
           (p->frozen.throat.temperature *
            p->frozen.exit.pressure *
            p->frozen.exit.velocity ));
    printf("C* (m/s)   : %8.3f \t %8.3f\n", 
           p->frozen.chamber.pressure * p->frozen.throat.aera_dotm,
           p->frozen.chamber.pressure * p->frozen.throat.aera_dotm);
    printf("Cf         : %8.3f \t %8.3f\n",
           p->frozen.throat.velocity /
           (p->frozen.chamber.pressure * p->frozen.throat.aera_dotm),
           p->frozen.exit.velocity /
           (p->frozen.chamber.pressure * p->frozen.throat.aera_dotm));
    printf("Ivac (m/s) : %8.3f \t %8.3f\n",
           p->frozen.throat.velocity + p->frozen.throat.pressure
           * p->frozen.throat.aera_dotm,
           p->frozen.exit.velocity + p->frozen.exit.pressure
           * p->frozen.exit.aera_dotm);
    printf("Isp (m/s)  : %8.3f \t %8.3f\n",
           p->frozen.throat.velocity, p->frozen.exit.velocity);


  }

  if (p->equilibrium_ok)
  {
    fprintf(outputfile,
            "\n--- Shifting equilibrium performance characteristics. ---\n");

    
    printf("Ae/At      : %8.3f \t %8.3f\n", 1.0,
           (p->equilibrium.exit.temperature *
            p->equilibrium.throat.pressure *
            p->equilibrium.throat.velocity) /
           (p->equilibrium.throat.temperature *
            p->equilibrium.exit.pressure *
            p->equilibrium.exit.velocity ));
    printf("C* (m/s)   : %8.3f \t %8.3f\n", 
           p->equilibrium.chamber.pressure * p->equilibrium.throat.aera_dotm,
           p->equilibrium.chamber.pressure * p->equilibrium.throat.aera_dotm);
    printf("Cf         : %8.3f \t %8.3f\n",
           p->equilibrium.throat.velocity /
           (p->equilibrium.chamber.pressure * p->equilibrium.throat.aera_dotm),
           p->equilibrium.exit.velocity /
           (p->equilibrium.chamber.pressure * p->equilibrium.throat.aera_dotm));
    printf("Ivac (m/s) : %8.3f \t %8.3f\n",
           p->equilibrium.throat.velocity + p->equilibrium.throat.pressure
           * p->equilibrium.throat.aera_dotm,
           p->equilibrium.exit.velocity + p->equilibrium.exit.pressure
           * p->equilibrium.exit.aera_dotm);
    printf("Isp (m/s)  : %8.3f \t %8.3f\n",
           p->equilibrium.throat.velocity, p->equilibrium.exit.velocity);
    

  }

  return 0;
}
*/





