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

#ifndef num_h
#define num_h


/**************************************************************
FUNCTION: This is a function to solve a matrix of linear
          equation. It use the Gauus elimination method as
	  explain in 'Advanced engineering mathematics' bye
	  Erwin Kreyszig.

PARAMETER: **matrix is a pointer to an array of [neq]x[neq+1]
           which is an augmented matrix.
	   *solution is an array that will contain the solution
	   if there is one and junk if there is no solution. It
	   must be allocated.
	   neq is the number of unknoen in the system

OUTPUT: 0 on success, -1 on failure

COMMENTS: the number of operation can be estimated by 2*n^3/3

AUTHOR:  Antoine Lefebvre

DATE: February 6, 2000
***************************************************************/
int gauss(double **matrix, double *solution, int neq);


/**************************************************************
FUNCTION: This function solve system of linear equation with
          the LU-Factorisation method as explain in
	  'Advanced engineering mathematics' bye Erwin Kreyszig.

THEORY:   To solve the system Ax=b , we could decompose A such
          that A = LU. L is a lower triangular matrix and
	  U is an upper triangular matrix. We could then solve
	  the system bye substitution cause Ly=b and Ux=y.

PARAMETER: Same as for gauss

OUTPUT: 0 on success, -1 on failure

COMMENTS: the number of operation is estimate by n^3/3, about
          half as many as the Gauss elimination. 
	  Should be interesting to compare the speed in case
	  of big matrix

AUTHOR:  Antoine Lefebvre

DATE: February 6, 2000
***************************************************************/
int lu(double **matrix, double *solution, int neq);


/**************************************************************
FUNCTION: This function print the coefficient of the matrix to
          the screen. 

PARAMETER: Same as for gauss
***************************************************************/
int print_matrix(double **matrix, int neq);

/* same thing but for square matrix instead of (N)x(N-1) */
int print_square_matrix(double **matrix, int neq);

/**************************************************************
FUNCTION: This function print the contents of the vector to
          the screen. 

PARAMETER: Same as for gauss
***************************************************************/
int print_vec(double *vec, int neq);

#endif
