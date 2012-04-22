#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#ifdef _MSC_VER
#define INLINE __inline
#else
#define INLINE inline
#endif

#ifdef _MSC_VER

INLINE static double fmax(double a, double b)
{
	return a>=b?a:b;
}

INLINE static double fmin(double a, double b)
{
	return a<=b?a:b;
}

#endif

#endif
