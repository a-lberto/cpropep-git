#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef GCC
#include <unistd.h>
#else
#include "getopt.h"
#endif

#include "thermo.h"
#include "load.h"

#define MAX_LINE 500
#define ARG_LENGTH 32

#define PROMPT "PROP>"
#define VERSION "1.0"

/* Dimensionless Enthalpy */
#define ENTHALPY(n, T) (enthalpy_0(n, T) - ((thermo_list + n)->heat -\
                       (thermo_list + n)->dho)/(R*T))

/* Dimensionless Entropy */
#define ENTROPY(n, T)  entropy_0(n, T)

/* Dimensionless Internal Energy */
#define INTERNAL_ENERGY(n, T) (enthalpy_0(n, T) - ((thermo_list + n)->heat -\
                               (thermo_list + n)->dho)/(R*T) - 1)

/* Dimensionless Specific Heat */
#define SPECIFIC_HEAT(n, T) specific_heat_0(n, T)

/* Dimensionless gibbs function */
#define GIBBS(n, T) gibbs_0(n, T)

#define MOLAR_MASS(n) (thermo_list + n)->weight

#define NAME(n) (thermo_list + n)->name

int global_verbose = 0;


typedef struct _command_
{
  char name[64];
  int  num_arg;
  void (*function)(char **arg);
  char description[512];
  char syntax[256];
} command_t;


void help(char **arg);
void molarMass(char **arg);
void Entropy(char **arg);
void Gibbs(char **arg);
void Enthalpy(char **arg);
void Specific_heat(char **arg);
void Internal_energy(char **arg);
void Properties(char **arg);
void List(char **arg);
void Info(char **arg);

int converti(char *string, unsigned int *number);
int convertd(char *string, double *number);

command_t command_list[] = {
  {"properties",      2, (void*)Properties,
   "All properties.\n", "<molecule> <temp>"},
  {"enthalpy",        2, (void*)Enthalpy,
   "Enthalpy.\n", "<molecule> <temp>"},
  {"internal_energy", 2, (void*)Internal_energy,
   "Internal energy.\n", "<molecule> <temp>"},
  {"entropy",         2, (void*)Entropy,
   "Entropy.\n", "<molucule> <temp>"},
  {"specific_heat",   2, (void*)Specific_heat,
   "Specific heat (Cp).\n", "<molecule> <temp>"},
  {"gibbs",           2, (void*)Gibbs,
   "Gibbs free energy.\n", "<molecule> <temp>"},
  {"molar_mass",      1, (void*)molarMass,
   "Molar mass.\n", "<molucule>"},
  {"info",            1, (void*)Info,
   "Reference information.\n", "<molecule>"},
  {"list",            1, (void*)List,
   "List molecule code beginning by formula.\n", "<formula>"},
  {"quit",            0, (void*)exit, "Exit the program.\n", ""},
  {"bye",             0, (void*)exit, "Exit the program.\n", ""},
  {"exit",            0, (void*)exit, "Exit the program.\n", ""},
  {"help",            0, (void*)help, "Display help message.\n", ""},
  {"\0",              0,  NULL, "\0"}
};

void help(char **arg)
{
  int i = 0;
  printf("\n");
  printf("%-36s %s\n", "Command", "Description");
  printf("%-36s %s\n", "-------", "-----------");
  
  while (command_list[i].name[0] != '\0')
  {
    printf("%-15s %-20s %s", command_list[i].name, command_list[i].syntax,
           command_list[i].description);
    i++;
  }
  printf("\n");

  printf("Arguments\n");
  printf("---------\n");
  printf("<molecule>      Number representing a specific molecule in\n");
  printf("                the database.\n");
  printf("<temp>          Temperature in deg K.\n");
  printf("<formula>       Chemical notation (ex: Fe2Cl)\n");
  printf("\n");
}

void usage(void)
{
  printf("Usage:\n");
  printf("        -d data_file\n");
  printf("        -fh\n\n");
  printf("Arguments:\n");
  printf("-d data_file         Path of the thermodynamics data file\n");
  printf("-f                   Supress introduction message\n");
  printf("-h                   Help message\n");
}

void List(char **arg)
{
  int i;

#ifdef BORLAND
  int j = 0;
#endif
  
  char formula[128];
  strncpy(formula, arg[0], 128);
  for (i = 0; i <  num_thermo; i++)
  {
    if ((STRNCASECMP(NAME(i), formula, strlen(formula)-1)) == 0)
    {
#ifdef BORLAND
      j++;
      if (j%22 == 0)
      {
        printf("(Press any key to continue)");
        getchar();
      }
#endif
      printf("%4d: %s\n", i, NAME(i));
      
    }
  }
  return;
  
}

void Info(char **arg)
{
  unsigned int n;
  /* molecule temperature pressure */
  if (converti(arg[0], &n) != 0)
  {
    printf("Error converting number.\n");
    return;
  }
  if (n > num_thermo)
  {
    printf("Value out of range.\n");
    return;
  }
  printf("\nInformation about %s data.\n", NAME(n));
  printf("---------------------------\n");
  printf("Reference: %s\n", (thermo_list + n)->comments);
  printf("Id       : %s\n\n", (thermo_list + n)->id);
  return;
}

void Properties(char **arg)
{
  unsigned int n;
  double T;
  
  /* molecule temperature pressure */
  if (converti(arg[0], &n) != 0)
  {
    printf("Error converting number.\n");
    return;
  }
  if (n > num_thermo)
  {
    printf("Value out of range.\n");
    return;
  }

  if (convertd(arg[1], &T) != 0)
  {
    printf("Error converting number.\n");
    return;
  }

  printf("\nProperties of %s at %.2f deg K\n", NAME(n), T);
  printf("---------------------------------\n");
  printf("    Enthalpy      -  Internal energy     -   Entropy\n");
  printf("%9.3f %-11s %9.3f %-9s %9.3f %-10s\n",
         ENTHALPY(n, T), "(-)",
         INTERNAL_ENERGY(n, T), "(-)",
         ENTROPY(n, T), "(-)");
  printf("%9.3f %-11s %9.3f %-9s %9.3f %-10s\n",
         ENTHALPY(n, T)*R*T/1000, "(kJ/mol)", 
         INTERNAL_ENERGY(n, T)*R*T/1000, "(kJ/mol)",
         ENTROPY(n, T)*R, "(kJ/kmol-K)");
  printf("%9.3f %-11s %9.3f %-9s %9.3f %-10s\n\n",
         ENTHALPY(n, T)*R*T/MOLAR_MASS(n), "(kJ/kg)",
         INTERNAL_ENERGY(n, T)*R*T/MOLAR_MASS(n), "(kJ/kg)",
         ENTROPY(n, T)*R/MOLAR_MASS(n), "(kJ/kg-K)");

  printf("  Specific Heat   -  Gibbs free energy\n");
  printf("%9.3f %-11s %9.3f %-9s\n",
         SPECIFIC_HEAT(n, T), "(-)",
         GIBBS(n, T), "(-)");
  printf("%9.3f %-11s %9.3f %-9s\n",
         SPECIFIC_HEAT(n, T)*R, "(kJ/kmol-K)",
         GIBBS(n, T)*R*T/1000, "(kJ/mol)");
  printf("%9.3f %-11s %9.3f %-9s\n\n",
         SPECIFIC_HEAT(n, T)*R/MOLAR_MASS(n), "(kJ/kg-K)",
         GIBBS(n, T)*R*T/MOLAR_MASS(n), "(kJ/kg)");

}

void Internal_energy(char **arg)
{
  unsigned int n;
  double T;
  
  /* molecule temperature pressure */
  if (converti(arg[0], &n) != 0)
  {
    printf("Error converting number.\n");
    return;
  }
  if (n > num_thermo)
  {
    printf("Value out of range.\n");
    return;
  }

  if (convertd(arg[1], &T) != 0)
  {
    printf("Error converting number.\n");
    return;
  }

  printf("\nInternal Energy of %s at %.2f\n", NAME(n), T);
  printf("%9.3f (-)\n%9.3f (kJ/mol)\n %9.3f (kJ/kg)\n\n",
         INTERNAL_ENERGY(n, T),
         INTERNAL_ENERGY(n, T)*R*T/1000,
         INTERNAL_ENERGY(n, T)*R*T/MOLAR_MASS(n));
  
  return;
 
}

void Entropy(char **arg)
{
  unsigned int n;
  double T;
  
  /* molecule temperature pressure */
  if (converti(arg[0], &n) != 0)
  {
    printf("Error converting number.\n");
    return;
  }
  if (n > num_thermo)
  {
    printf("Value out of range.\n");
    return;
  }

  if (convertd(arg[1], &T) != 0)
  {
    printf("Error converting number.\n");
    return;
  }

  printf("\nEntropy of %s at %.2f deg K\n", NAME(n), T);
  printf("%9.3f (-)\n%9.3f kJ/kmol-K\n%9.3f kJ/kg-K\n\n",
         ENTROPY(n, T),
         ENTROPY(n, T)*R,
         ENTROPY(n, T)*R/MOLAR_MASS(n));
  
  return;
}

void Gibbs(char **arg)
{
  unsigned int n;
  double T;
  
  /* molecule temperature pressure */
  if (converti(arg[0], &n) != 0)
  {
    printf("Error converting number.\n");
    return;
  }
  if (n > num_thermo)
  {
    printf("Value out of range.\n");
    return;
  }

  if (convertd(arg[1], &T) != 0)
  {
    printf("Error converting number.\n");
    return;
  }

  printf("\nGibbs function of %s at %.2f deg K\n", NAME(n), T);
  printf("%9.3f (-)\n%9.3f (kJ/mol)\n%9.3f (kJ/kg)\n\n",
         GIBBS(n, T),
         GIBBS(n, T)*R*T/1000,
         GIBBS(n, T)*R*T/MOLAR_MASS(n));
  
  return;
}

void Enthalpy(char **arg)
{
  unsigned int n;
  double T;
  
  /* molecule temperature pressure */
  if (converti(arg[0], &n) != 0)
  {
    printf("Error converting number.\n");
    return;
  }
  if (n > num_thermo)
  {
    printf("Value out of range.\n");
    return;
  }

  if (convertd(arg[1], &T) != 0)
  {
    printf("Error converting number.\n");
    return;
  }

  printf("\nEnthalpy of %s at %.2f deg K\n", NAME(n), T);
  printf("%9.3f (-)\n%9.3f (kJ/mol)\n%9.3f (kJ/kg)\n\n",
         ENTHALPY(n, T),
         ENTHALPY(n, T)*R*T/1000,
         ENTHALPY(n, T)*R*T/MOLAR_MASS(n));
  
  return;
}

void Specific_heat(char **arg)
{
  unsigned int n;
  double T;
  
  /* molecule temperature pressure */
  if (converti(arg[0], &n) != 0)
  {
    printf("Error converting number.\n");
    return;
  }
  if (n > num_thermo)
  {
    printf("Value out of range.\n");
    return;
  }

  if (convertd(arg[1], &T) != 0)
  {
    printf("Error converting number.\n");
    return;
  }

  printf("\nSpecific heat (Cp) of %s at %.2f\n", NAME(n), T);
  printf("%9.3f (-)\n%9.3f (kJ/kmol-K)\n%9.3f (kJ/kg-K)\n\n",
         SPECIFIC_HEAT(n, T),
         SPECIFIC_HEAT(n, T)*R,
         SPECIFIC_HEAT(n, T)*R/MOLAR_MASS(n));
  
  return;
}

void molarMass(char **arg)
{
  unsigned int n;

  if (converti(arg[0], &n) != 0)
  {
    printf("Error converting number.\n");
    return;
  }
  
  if (n > num_thermo)
  {
    printf("Value out of range.\n");
    return;
  }

  printf("\nMolar mass of %s\n", NAME(n));
  printf("%9.4f g/mol\n\n", MOLAR_MASS(n));

  return;
}


void welcome(void)
{
  printf("PROP, version %s (GPL License).\n", VERSION);
  printf("Ideal Gas Thermochemical Properties.\n");
  printf("Interface to the NASA Gleen Research Center Database\n");
  printf("provide by Bonnie J. McBride.\n");
  printf("Copyright (C) 2000 Antoine Lefebvre.\n");
  printf("For command list, type `help'.\n\n");
}

int converti(char *string, unsigned int *number)
{

  char *ptr;
  
  if (string[0] == '\0')
  {
    return -1;
  }

  errno = 0; 
  *number = strtoul(string, &ptr, 10);
  
  if (errno == ERANGE)
  {
    return -1;
  }

  if (ptr[0] == string[0])
  {
    return -1;
  }

  return 0;
}

int convertd(char *string, double *number)
{
  char *ptr;
  
  if (string[0] == '\0')
  {
    return -1;
  }

  errno = 0; 
  *number = strtod(string, &ptr);
  
  if (errno == ERANGE)
  {
    return -1;
  }

  if (ptr[0] == string[0])
  {
    return -1;
  }

  return 0;
}


int main(int argc, char *argv[])
{

  int i, j, n, c;
  bool found;
  bool ok;
  bool msg = true;
  
  char line[MAX_LINE];
  char *ptr;
  char *running;
  char *token;

  char **arg;
  
  char thermo_file[FILENAME_MAX] = "thermo.dat";

  /* parse the command line */
  while (1)
  {
    c = getopt(argc, argv, "fhd:");

    if (c == EOF)
      break;
    
    switch (c)
    {
      case 'd':
          strncpy(thermo_file, optarg, FILENAME_MAX);
          break;
      case 'f':
          msg = false;
          break;
      case 'h':
          usage();
          return 0;
    }
  }
  

  if (msg)
    welcome();
  
  if (load_thermo (thermo_file) < 0)
  {
    printf("Error loading thermo data file: %s\n", thermo_file);
    return -1;
  }
  
  while (1)
  {
    printf("%s", PROMPT);
    
    if (fgets(line, MAX_LINE, stdin) == NULL)
    {
      printf("\n");
      return 0;
    }

    if (line[0] != '\n')
    {
      i = 0;
      found = false;

      /* duplicate the string */
      if ((running = strdup(line)) == NULL)
      {
        printf("Not enough memory.\n");
        return -1;
      }
      /* keep the adress in memory to free it later */
      ptr = running;
      /* first get the command name */
      token = strtok(running, " \n");
      
      while ((command_list[i].name[0] != '\0') && (found == false))
      {
        
        if (STRNCASECMP(command_list[i].name, token,
                        strlen(token)) == 0)
        {
          found = true;
          ok = true;
          
          /* number of arguments for this command */
          n = command_list[i].num_arg;
          
          /* allocate the memory for the arguments */
          arg = (char **) malloc(n * sizeof(char*));
          for (j = 0; j < n; j++)
            arg[j] = (char *) malloc(ARG_LENGTH * sizeof(char));
          
          /* get each arguments */
          for (j = 0; j < n; j++)
          {
            if ((token = strtok(NULL, " ")) == NULL)
            {
              printf("Not enough arguments.\n");
              ok = false;
              break;
            }
            /* copy it */
            if (ok)
              strncpy(arg[j], token, ARG_LENGTH);
          }
          
          if (ok)
          {
            /* verify if there is more arguments than needed */
            if ((token = strtok(NULL, " ")) != NULL)
              printf("Too much arguments. Need only %d.\n", n);
            
            /* execute the command */
            command_list[i].function(arg);
          }
          
          /* free memory */
          for (j = 0; j < n; j++)
            free(arg[j]);
          free(arg);
        }
        i++;
      }
      free(ptr);
    
      if (found == false)
      {
        printf("?Invalid command\n");
      }
    }
  }
  return 0;
}
