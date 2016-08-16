#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#ifndef _WIN32
#define TEXT(X) X
typedef char TCHAR;
typedef TCHAR* LPTSTR;
#else
#include <windows.h>
#endif

#ifdef _MSC_VER
#define INLINE __inline
#else
#define INLINE inline
#endif

#ifdef CRT_NO_FMAX
/* msvcrt does not have fmin() and fmax() */

INLINE static double fmax(double a, double b)
{
	return a>=b?a:b;
}

INLINE static double fmin(double a, double b)
{
	return a<=b?a:b;
}
#endif


extern void (*mom_trace)(TCHAR *);
extern void (*mom_error)(TCHAR *);

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define C0	3e8
#define Z0	(120 * M_PI)
#define MU0	(Z0 / C0)
#define EPS0	(1 / (Z0*C0))

#endif
