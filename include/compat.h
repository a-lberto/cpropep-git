#ifndef COMPAT_H
#define COMPAT_H

/*
  Checking for _MSC_VER will detect whether MSVC++ is being used.
  I don't know of the other compiler flags, so others will want to
  add to this for their own compilers.
  
  Mark Pinese 24/4/2000
*/



#ifdef _MSC_VER
  /* MSVC++ 6.0 Std */
#define STRNCASECMP		_strnicmp
	
#ifdef _DEBUG
#include <crtdbg.h>
#endif /* defined (_DEBUG) */

#ifndef __cplusplus
typedef enum
{
  false = 0,
  true = 1
} bool;
#endif /* !defined (__cplusplus) */

#endif /* define _MSC_VER */

#ifdef GCC

#define STRNCASECMP		strncasecmp
#define __min(a, b) ( (a) <= (b) ? (a) : (b))
#define __max(a, b) ( (a) >= (b) ? (a) : (b))

typedef enum
{
  false = 0,
  true  = 1
} bool;

#endif /* define GCC */

#ifdef BORLAND

int StrNCaseCmp(const char *s1, const char *s2, size_t sz);

#define STRNCASECMP StrNCaseCmp
#define __min(a, b) ( (a) <= (b) ? (a) : (b))
#define __max(a, b) ( (a) >= (b) ? (a) : (b))

typedef enum
{
  false = 0,
  true  = 1
} bool;

#endif /* define BORLAND */
    

#endif	/* !defined(COMPAT_H) */
