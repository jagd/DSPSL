#ifndef __MOM_H__
#define __MOM_H__

#include "mom_mesh.h"
#include "md.h"

double mom(
	struct MeshConfig *conf,
	/*
	   charge density vector for freespace
	   with cdfs[0] := only free charge
	        cdfs[1] := all charge
	   index range upto conf->index[ID_STRIP_END];
	*/
	struct MD *cdfs[], /* default NULL */
	/*
	   charge density vector with dielectric
	   with cdd[0] := only free charge
	        cdd[1] := all charge
	   index range upto conf->index[ID_MESH_CELLS];
	*/
	struct MD *cdd[], /* default NULL */
	double capacity[2] /* default NULL */
	);

#endif
