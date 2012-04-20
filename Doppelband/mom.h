#ifndef __MOM_H__
#define __MOM_H__

#include "mom_mesh.h"
#include "md.h"

struct MD* mom_matrix_new(struct MeshConfig *conf);
struct MD* calc_charge(
		struct MeshConfig *conf,
		struct MD* a,
		double charge[/* 2 */] );

#endif
