#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mom_mesh.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define C0	3e8
#define Z0	(120 * M_PI)
#define MU0	(Z0 / C0)
#define EPS0	(1 / (Z0*C0))


