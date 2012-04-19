#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mom_mesh.h"

#define error printf
#define log printf
#define INLINE inline

struct Line_1D {
	double left, right;
};


struct MeshState  { /* local state */
	double w[2], h, offset, w_port_ext;
	struct Line_1D port, strip[2];
	struct MeshConfig *conf;
	double mesh_step;
	int n; /* current cell index */
};


void mesh_free(struct MeshConfig *m)
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
	double width;

	/* calculating the coordinates */
	width = fmax(s->w[0]*0.5, s->w[1]*0.5 - s->offset) + s->offset
				+ fmax(s->w[0]*0.5 - s->offset, s->w[1]*0.5);

	s->port.left = -width*0.5 - s->w_port_ext;
	s->port.right = width*0.5 + s->w_port_ext;

	s->strip[0].left = width*0.5
			- fmax(s->w[0]*0.5, s->w[1]*0.5 - s->offset)
			- s->w[0]*0.5;
	s->strip[0].right = width*0.5
			- fmax(s->w[0]*0.5, s->w[1]*0.5 - s->offset)
			+ s->w[0]*0.5;

	s->strip[1].left = -width*0.5
			+ fmax(s->w[0]*0.5 - s->offset, s->w[1]*0.5)
			- s->w[1]*0.5;
	s->strip[1].right = -width*0.5
			+ fmax(s->w[0]*0.5 - s->offset, s->w[1]*0.5)
			+ s->w[1]*0.5;

}


static INLINE void mesh_auto_predict(struct MeshState *s)
{
	int max_cells;
	/* define the general mesh length */
	s->mesh_step = fmin(s->h*0.2, (s->w[0] + s->w[1])*0.05);

	/* consider the edge refine */
	/*
	   12 edges, each has 2 extra mesh cell ==> 24 extra cells
	*/
	max_cells = ceil(
		2 * (s->port.right - s->port.left)/s->mesh_step + 24 );

	s->conf = (struct MeshConfig*)malloc(sizeof(struct MeshConfig));
	s->conf->mesh = (struct Cell_1D*)malloc(
			sizeof(struct Cell_1D) * max_cells);

#ifdef MOM_MESH_ENABLE_DEBUG
	log("Predicted number of cells = %d\n", max_cells);
#endif

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
	s->conf->index[ID_STRIP0_START] = s->n;

	segs = floor(s->w[0] / s->mesh_step);
	if (segs > 0) { /* reserve some space for edge refining */
		--segs;
	}

	rest = s->w[0] - (segs*s->mesh_step);

	l_short = rest * (3.0/16.0);
	l_long = rest * (5.0/16.0);

	x = s->strip[0].left;

	/* refining left edge*/
	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	/* the normal segments */
	cache = s->mesh_step * 0.5;
	for (i = 0; i < segs; ++i) {
		x += cache;
		s->conf->mesh[s->n].centre = x;
		s->conf->mesh[s->n].hw = cache;
		x += cache;
		s->n++;
	}

	/* refining right edge*/
	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;


	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
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

	/* meshing the STRIP1 */
	s->conf->index[ID_STRIP1_START] = s->n;

	segs = floor(s->w[1] / s->mesh_step);
	if (segs > 0) { /* reserve some space for edge refining */
		--segs;
	}

	rest = s->w[1] - (segs*s->mesh_step);

	l_short = rest * (3.0/16.0);
	l_long = rest * (5.0/16.0);

	x = s->strip[1].left;

	/* refining left edge*/
	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	/* the normal segments */
	cache = s->mesh_step * 0.5;
	for (i = 0; i < segs; ++i) {
		x += cache;
		s->conf->mesh[s->n].centre = x;
		s->conf->mesh[s->n].hw = cache;
		x += cache;
		s->n++;
	}

	/* refining right edge*/
	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;
}


static INLINE void mesh_dielectric0_left(struct MeshState *s)
{
	int i;

	int segs; /* segments in general mesh step */
	double x;
	double rest; /* the total length of refined edge for each element */
	double l_short, l_long;
	double cache;
	double w;

	s->conf->index[ID_DIELECTRIC0_START] = s->n;

	w = s->strip[0].left - s->port.left;

	segs = floor(w / s->mesh_step);
	if (segs > 0) { /* reserve some space for edge refining */
		--segs;
	}

	rest = w - (segs*s->mesh_step);

	l_short = rest * (3.0/16.0);
	l_long = rest * (5.0/16.0);

	x = s->port.left;

	/* refining left edge*/
	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	/* the normal segments */
	cache = s->mesh_step * 0.5;
	for (i = 0; i < segs; ++i) {
		x += cache;
		s->conf->mesh[s->n].centre = x;
		s->conf->mesh[s->n].hw = cache;
		x += cache;
		s->n++;
	}

	/* refining right edge*/
	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;
}


static INLINE void mesh_dielectric0_right(struct MeshState *s)
{
	int i;

	int segs; /* segments in general mesh step */
	double x;
	double rest; /* the total length of refined edge for each element */
	double l_short, l_long;
	double cache;
	double w;

	w = s->port.right - s->strip[0].right;

	segs = floor(w / s->mesh_step);
	if (segs > 0) { /* reserve some space for edge refining */
		--segs;
	}

	rest = w - (segs*s->mesh_step);

	l_short = rest * (3.0/16.0);
	l_long = rest * (5.0/16.0);

	x = s->strip[0].right;

	/* refining left edge*/
	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	/* the normal segments */
	cache = s->mesh_step * 0.5;
	for (i = 0; i < segs; ++i) {
		x += cache;
		s->conf->mesh[s->n].centre = x;
		s->conf->mesh[s->n].hw = cache;
		x += cache;
		s->n++;
	}

	/* refining right edge*/
	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;
}


static INLINE void mesh_dielectric1_left(struct MeshState *s)
{
	int i;

	int segs; /* segments in general mesh step */
	double x;
	double rest; /* the total length of refined edge for each element */
	double l_short, l_long;
	double cache;
	double w;

	s->conf->index[ID_DIELECTRIC1_START] = s->n;

	w = s->strip[1].left - s->port.left;

	segs = floor(w / s->mesh_step);
	if (segs > 0) { /* reserve some space for edge refining */
		--segs;
	}

	rest = w - (segs*s->mesh_step);

	l_short = rest * (3.0/16.0);
	l_long = rest * (5.0/16.0);

	x = s->port.left;

	/* refining left edge*/
	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	/* the normal segments */
	cache = s->mesh_step * 0.5;
	for (i = 0; i < segs; ++i) {
		x += cache;
		s->conf->mesh[s->n].centre = x;
		s->conf->mesh[s->n].hw = cache;
		x += cache;
		s->n++;
	}

	/* refining right edge*/
	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;
}


static INLINE void mesh_dielectric1_right(struct MeshState *s)
{
	int i;

	int segs; /* segments in general mesh step */
	double x;
	double rest; /* the total length of refined edge for each element */
	double l_short, l_long;
	double cache;
	double w;

	w = s->port.right - s->strip[1].right;

	segs = floor(w / s->mesh_step);
	if (segs > 0) { /* reserve some space for edge refining */
		--segs;
	}

	rest = w - (segs*s->mesh_step);

	l_short = rest * (3.0/16.0);
	l_long = rest * (5.0/16.0);

	x = s->strip[1].right;

	/* refining left edge*/
	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	/* the normal segments */
	cache = s->mesh_step * 0.5;
	for (i = 0; i < segs; ++i) {
		x += cache;
		s->conf->mesh[s->n].centre = x;
		s->conf->mesh[s->n].hw = cache;
		x += cache;
		s->n++;
	}

	/* refining right edge*/
	cache = (l_long * 0.5);
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;

	cache = l_short * 0.5;
	x += cache;
	s->conf->mesh[s->n].centre = x;
	s->conf->mesh[s->n].hw = cache;
	x += cache;
	s->n++;
}

/*
   near the edge, mesh will be finer,
   the general mesh_step will be divided into  3/8 and 5/8
*/
struct MeshConfig* mesh_new(
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

	mesh_dielectric0_left(&s);
	mesh_dielectric0_right(&s);
	mesh_dielectric1_left(&s);
	mesh_dielectric1_right(&s);

	s.conf->index[ID_MESH_CELLS] = s.n;

#ifdef MOM_MESH_ENABLE_DEBUG
	log("Port = [ %lf\t%lf ]\n", s.port.left, s.port.right);
	log("Strip[0] = [ %lf\t%lf ]\n", s.strip[0].left, s.strip[0].right);
	log("Strip[1] = [ %lf\t%lf ]\n", s.strip[1].left, s.strip[1].right);
	log("Used number of cells = %d\n", s.conf->index[ID_MESH_CELLS]);
	log("Number of cells for strip[0] = %d\n",
		s.conf->index[ID_STRIP0_END]
		- s.conf->index[ID_STRIP0_START]);
	log("Number of cells for strip[1] = %d\n",
		s.conf->index[ID_STRIP1_END]
		- s.conf->index[ID_STRIP1_START]);
	log("Number of cells for dielectric[0] = %d\n",
		s.conf->index[ID_DIELECTRIC0_END]
		- s.conf->index[ID_DIELECTRIC0_START]);
	log("Number of cells for dielectric[1] = %d\n",
		s.conf->index[ID_DIELECTRIC1_END]
		- s.conf->index[ID_DIELECTRIC1_START]);
#endif
	return s.conf;
}
