#ifndef RETURN_H
#define RETURN_H

/* Codes used in some functions by Antoine Lefebvre */
#define  SUCCESS  0
#define  ERROR   -1

/*
  Return codes
  Mark Pinese 24/4/2000
*/

#define ERR_BASE	    -100
#define ERR_MALLOC	  ERR_BASE - 1
#define ERR_FOPEN		  ERR_BASE - 2
#define ERR_EOF			  ERR_BASE - 3
#define ERR_NOT_ALLOC	ERR_BASE - 4

#endif	/* !defined(RETURN_H) */
