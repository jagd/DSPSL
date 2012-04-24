#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "mom_mesh.h"

#ifdef MOM_MESH_ENABLE_DEBUG
#define DEBUG_BUF_SIZE 1000
static TCHAR debug_buf[DEBUG_BUF_SIZE];
#endif

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


/* return: success 0 else non-zero */
static INLINE int mesh_auto_predict(struct MeshState *s, double h)
{
	int max_cells;
	/* define the general mesh length */
	s->mesh_step = fmin(h*0.15, (s->w[0] + s->w[1])*0.025);

	/* consider the edge refine */
	/*
	   12 edges, each has 2 extra mesh cell ==> 24 extra cells
	*/
	max_cells = (int)ceil(
		2 * (s->port.right - s->port.left)/s->mesh_step + 24 );

	s->conf = (struct MeshConfig*)malloc(sizeof(struct MeshConfig));
	s->conf->mesh = (struct Cell_1D*)malloc(
			sizeof(struct Cell_1D) * max_cells);
	if (s->conf == 0 || s->conf->mesh == 0) {
		mom_error(TEXT("Create configuration failed"));
		if (s->conf) {
			mesh_free(s->conf);
		}
		return 1;
	}

#ifdef MOM_MESH_ENABLE_DEBUG
	snprintf(debug_buf, DEBUG_BUF_SIZE
			,TEXT("Predicted number of cells = %d\n")
			, max_cells);
	mom_trace(debug_buf);
#endif
	return 0;
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

	segs = (int)floor(s->w[0] / s->mesh_step);
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

	segs = (int)floor(s->w[1] / s->mesh_step);
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

	segs = (int)floor(w / s->mesh_step);
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

	segs = (int)floor(w / s->mesh_step);
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

	segs = (int)floor(w / s->mesh_step);
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

	segs = (int)floor(w / s->mesh_step);
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
		double h, /* only for determination of the mesh length */
		double eps_r /* epsilon_r */
		)
{
	struct MeshState s;

	/* check inputed parameters */
	if (w0 <= 0 || w1 <= 0 || w_port_ext <= 0 || h <= 0) {
		mom_error(TEXT("generate_mesh(): ")
				TEXT("Illegal geometric length"));
		return NULL;
	}

	s.w[0] = w0;
	s.w[1] = w1;
	s.offset = offset;
	s.w_port_ext = w_port_ext;

	mesh_calc_coord(&s);
	if (mesh_auto_predict(&s, h) != 0) {
		mom_error(TEXT("generate_mesh(): ")
				TEXT("not enough memory"));
		return NULL;
	}

	s.n = 0;

	mesh_strip0(&s);
	mesh_strip1(&s);

	/* too small eps_r cause singularity */
	if (eps_r > EPS_MIN) {
		mesh_dielectric0_left(&s);
		mesh_dielectric0_right(&s);
		mesh_dielectric1_left(&s);
		mesh_dielectric1_right(&s);
	} else {
		if (eps_r != 1.0) {
			mom_trace(TEXT("WARNING: ")
				TEXT("epsilon should be either 1.0 ")
				TEXT("or much bigger than 1"));
			mom_trace(TEXT("Calculation with epsilon_r = 1.0"));
			eps_r = 1.0;
		}


		s.conf->index[ID_DIELECTRIC0_START] = s.n;
		s.conf->index[ID_DIELECTRIC0_END] = s.n;
		s.conf->index[ID_DIELECTRIC1_START] = s.n;
		s.conf->index[ID_DIELECTRIC1_END] = s.n;
	}

	s.conf->index[ID_MESH_CELLS] = s.n;
	s.conf->h = h;
	s.conf->eps_r = eps_r;

#ifdef MOM_MESH_ENABLE_DEBUG
	snprintf(debug_buf, DEBUG_BUF_SIZE
			, TEXT("Used number of cells = %d")
			, s.conf->index[ID_MESH_CELLS]);
	mom_trace(debug_buf);

	snprintf(debug_buf, DEBUG_BUF_SIZE
			, TEXT("Number of cells for strip[0] = %d")
			, s.conf->index[ID_STRIP0_END]
				- s.conf->index[ID_STRIP0_START]);
	mom_trace(debug_buf);

	snprintf(debug_buf, DEBUG_BUF_SIZE
			, TEXT("Number of cells for strip[1] = %d")
			, s.conf->index[ID_STRIP1_END]
				- s.conf->index[ID_STRIP1_START]);
	mom_trace(debug_buf);

	snprintf(debug_buf, DEBUG_BUF_SIZE
			, TEXT("Number of cells for dielectric[0] = %d")
	      		, s.conf->index[ID_DIELECTRIC0_END]
				- s.conf->index[ID_DIELECTRIC0_START]);
	mom_trace(debug_buf);

	snprintf(debug_buf, DEBUG_BUF_SIZE
			, TEXT("Number of cells for dielectric[1] = %d")
			, s.conf->index[ID_DIELECTRIC1_END]
				- s.conf->index[ID_DIELECTRIC1_START]);
	mom_trace(debug_buf);

	snprintf(debug_buf, DEBUG_BUF_SIZE
			, TEXT("    Port = [ %le\t%le ]")
			, s.port.left, s.port.right);
	mom_trace(debug_buf);

	snprintf(debug_buf, DEBUG_BUF_SIZE
			, TEXT("Strip[0] = [ %le\t%le ]")
			, s.strip[0].left, s.strip[0].right);
	mom_trace(debug_buf);

	snprintf(debug_buf, DEBUG_BUF_SIZE
			,TEXT("Strip[1] = [ %le\t%le ]")
			, s.strip[1].left, s.strip[1].right);
	mom_trace(debug_buf);
#endif

	return s.conf;
}
