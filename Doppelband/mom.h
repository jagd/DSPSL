#ifndef __MOM_H__
#define __MOM_H__

#include "mom_mesh.h"
#include "md.h"

struct MD* calc_charge(
		struct MeshConfig *conf,
		double charge[/* 2 */] );

#endif
