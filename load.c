#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "equilibrium.h"
#include "load.h"


#define CAL_TO_JOULE 4.1868
#define CLBS_TO_CCM  27.679905

extern propellant_t	*propellant_list;
extern thermo_t	        *thermo_list;

extern char symb[][3];


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

  if (global_verbose)
  {
    printf("Loading thermo data file...");
    fflush(stdout);
  }

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
          }
        }
		    
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



        /* I'm not quite sure it is necessary */
        
        /* if the value is 0 and this is the same substance as
           the previous one but in a different state ... */
        if ( (thermo_list + i)->heat == 0.0 )
        {
          l = 1;
            
          for (j = 0; j < 5; j++)
          {
            if (!((thermo_list+i)->coef[j] == (thermo_list+i-1)->coef[j] &&
                  (thermo_list+i)->elem[j] == (thermo_list+i-1)->elem[j]))
              l = 0;
          }
          /* set to the same value as the previous one */
          if (l)
            (thermo_list+i)->heat = (thermo_list+i-1)->heat; 
        }
                    

            
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

  if (global_verbose)
    printf("%d species loaded.\n", i);
  
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

  if (global_verbose)
  {
    printf("Loading propellant data file...");
    fflush(stdout);
  }

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
      
      strncpy(tmp_ptr, buf_ptr + 69, 5);
      tmp[5] = '\0';		    
      propellant_list[i].heat = atoi(tmp) * CAL_TO_JOULE;
      
      strncpy(tmp_ptr, buf_ptr + 75, 5);
      tmp[5] = '\0';
      propellant_list[i].density = atof(tmp) * CLBS_TO_CCM;
      
      i++;
    }
  }  
  fclose(fd);     
  if (global_verbose)
    printf("%d species loaded.\n", i);
  return i;
}


void trim_spaces(char *str, unsigned int len)
{
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
