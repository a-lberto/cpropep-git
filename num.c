/* num.c  -  Library of numerical method                               */
/* Copyright (C) 2000                                                  */
/* Antoine Lefebvre <antoine.lefebvre@polymtl.ca                       */

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

//#include <time.h>

#include "libnum.h"


int main(void)
{

     /* for timing */
//     clock_t start, finish;


     const int neq = 3;

     int i = 0;

     double **mat;
     double *ans;

//     double mat[2][3];
//     double ans[2];

     mat = (double **)calloc(neq, sizeof(double*));
     for (i = 0; i < neq; i++)
     {
     mat[i] = (double *)calloc(neq+1, sizeof(double));
     }
     ans = (double *)calloc(neq, sizeof(double));

     mat[0][0] = 4;
     mat[0][1] = 10;
     mat[0][2] = -2;
     mat[0][3] = -20;
//     mat[0][4] = 54;

     mat[1][0] = -1;
     mat[1][1] = -15;
     mat[1][2] = 3;
     mat[1][3] = 30;
//     mat[1][4] = 0;

     mat[2][0] = 0;
     mat[2][1] = 25;
     mat[2][2] = -5;
     mat[2][3] = -50;
//     mat[2][4] = 2;

//     mat[3][0] = 43;
//     mat[3][1] = 0;
//     mat[3][2] = 76;
//     mat[3][3] = 0;
//     mat[3][4] = 77;

     print_matrix(mat, neq); 

     lu(mat, ans, neq);
  
//     gauss(mat, ans, neq);
     



     print_matrix(mat, neq);
     print_vec(ans, neq);

     //free(mat[1]);
     //free(mat[0]);
     free(mat);
     free(ans);
     return 0;

}

