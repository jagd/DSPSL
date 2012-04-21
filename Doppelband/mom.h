#ifndef __MOM_H__
#define __MOM_H__

#include "mom_mesh.h"
#include "md.h"

double mom(
	struct MeshConfig *conf,
	struct MD **cd_all, /* charge density for all charges */
	struct MD **cd_free);

#endif
