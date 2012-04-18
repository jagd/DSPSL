#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mom_mesh.h"

#define error printf
#define INLINE inline

struct Line_1D {
	double left, right;
};


struct MeshState  { /* local state */
	double w[2], h, offset, w_port_ext;
	double width; /* the total width */
	struct Line_1D port, strip[2];
	struct MeshConfig *conf;
	double mesh_step;
	int max_cells;
	int n; /* current cell index */
};


void MeshConfig_free(struct MeshConfig *m)
{
	if (m) {
		if (m->mesh) {
			free(m->mesh);
		}
		free(m);
	}
}


static INLINE void mesh_calc_coord(struct MeshState *s)
{
	/* calculating the coordinates */
	s->width = fmax(s->w[0]*0.5, s->w[1]*0.5 - s->offset) + s->offset
				+ fmax(s->w[0]*0.5 - s->offset, s->w[1]*0.5);

	s->port.left = -s->width*0.5 - s->w_port_ext;
	s->port.right = s->width*0.5 + s->w_port_ext;

	s->strip[0].left = s->width*0.5
			- fmax(s->w[0]*0.5, s->w[1]*0.5 - s->offset)
			- s->w[0]*0.5;
	s->strip[0].right = s->width*0.5
			- fmax(s->w[0]*0.5, s->w[1]*0.5 - s->offset)
			+ s->w[0]*0.5;

	s->strip[1].left = -s->width*0.5
			+ fmax(s->w[0]*0.5 - s->offset, s->w[1]*0.5)
			- s->w[1]*0.5;
	s->strip[1].right = -s->width*0.5
			+ fmax(s->w[0]*0.5 - s->offset, s->w[1]*0.5)
			+ s->w[1]*0.5;

}


static INLINE void mesh_auto_predict(struct MeshState *s)
{
	/* define the general mesh length */
	s->mesh_step = fmin(s->h*0.2, (s->w[0] + s->w[1])*0.05);

	/* consider the edge refine */
	/*
	   12 edges, each has 2 extra mesh cell ==> 24 extra cells
	*/
	s->max_cells = ceil(2*s->width/s->mesh_step + 24);

	s->conf = (struct MeshConfig*)malloc(sizeof(struct MeshConfig));
	s->conf->mesh = (struct Cell_1D*)malloc(
			sizeof(struct Cell_1D) * s->max_cells);

}


static INLINE void mesh_strip0(struct MeshState *s)
{
	int i;

	int segs; /* segments in general mesh step */
	double x;
	double rest; /* the total length of refined edge for each element */
	double l_short, l_long;
	double cache;

	/* meshing the STRIP0 */
	s->conf->index[ID_STRIP0_START] = 0;

	segs = floor(s->w[0] / s->mesh_step);
	if (segs > 0) { /* reserve some space for edge refining */
		--segs;
	}

	rest = s->w[0] - (segs*s->mesh_step);

	l_short = rest * (3/16);
	l_long = rest * (5/16);
	
	x = s->strip[0].left;

	/* refining left edge*/
	cache = l_short * 0.5;
	s->conf->mesh[s->n].centre = x + cache;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = (l_long * 0.5);
	s->conf->mesh[s->n].centre = x + cache;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	/* the normal segments */
	cache = s->mesh_step * 0.5;
	for (i = 0; i < segs; ++i) {
		s->conf->mesh[s->n].centre = x + cache;
		s->conf->mesh[s->n].hw = cache;
		x += cache;
		s->n++;
	}

	/* refining right edge*/
	cache = (l_long * 0.5);
	s->conf->mesh[s->n].centre = x + cache;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = l_short * 0.5;
	s->conf->mesh[s->n].centre = x + cache;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;
}

static INLINE void mesh_strip1(struct MeshState *s)
{
	int i;

	int segs; /* segments in general mesh step */
	double x;
	double rest; /* the total length of refined edge for each element */
	double l_short, l_long;
	double cache;

	/* meshing the STRIP0 */
	s->conf->index[ID_STRIP1_START] = 0;

	segs = floor(s->w[1] / s->mesh_step);
	if (segs > 0) { /* reserve some space for edge refining */
		--segs;
	}

	rest = s->w[1] - (segs*s->mesh_step);

	l_short = rest * (3/16);
	l_long = rest * (5/16);
	
	x = s->strip[1].left;

	/* refining left edge*/
	cache = l_short * 0.5;
	s->conf->mesh[s->n].centre = x + cache;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = (l_long * 0.5);
	s->conf->mesh[s->n].centre = x + cache;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	/* the normal segments */
	cache = s->mesh_step * 0.5;
	for (i = 0; i < segs; ++i) {
		s->conf->mesh[s->n].centre = x + cache;
		s->conf->mesh[s->n].hw = cache;
		x += cache;
		s->n++;
	}

	/* refining right edge*/
	cache = (l_long * 0.5);
	s->conf->mesh[s->n].centre = x + cache;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = l_short * 0.5;
	s->conf->mesh[s->n].centre = x + cache;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;
}


/*
   near the edge, mesh will be finer,
   the general mesh_step will be divided into  3/8 and 5/8
*/
struct MeshConfig* generate_mesh(
		double w0, /* strip0.width */
		double w1, /* strip1.width */
		double offset, /* offset of the centre of both strip */
		double w_port_ext, /* extend width in both side */
		double h /* only for determination of the mesh length */
		)
{
	struct MeshState s;

	/* check inputed parameters */
	if (w0 <= 0 || w1 <= 0 || w_port_ext <= 0 || h <= 0) {
		error("generate_mesh(): Illegal geometric length\n");
		return NULL;
	}

	s.w[0] = w0;
	s.w[1] = w1;
	s.h = h;
	s.offset = offset;
	s.w_port_ext = w_port_ext;

	mesh_calc_coord(&s);
	mesh_auto_predict(&s);

	s.n = 0;

	mesh_strip0(&s);
	mesh_strip1(&s);

	return NULL;
}
