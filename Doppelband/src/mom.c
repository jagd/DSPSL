#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mom_mesh.h"
#include "md.h"

void (*mom_trace)(char *) = NULL;

#define INLINE inline

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define C0	3e8
#define Z0	(120 * M_PI)
#define MU0	(Z0 / C0)
#define EPS0	(1 / (Z0*C0))

#define CONST_INV_2_PI_EPS0	(1.0 / (2.0*M_PI*EPS0))

#define INTEGRAL_STEPS 100

static INLINE void potential_strip0(
		struct MeshConfig *conf,
		struct MD *a,
		struct MD *k)
{
	double h = conf->h;
	int i;

	/* for every segments in strip[0]*/
	for (i = conf->index[ID_STRIP0_START];
		i < conf->index[ID_STRIP0_END]; ++i) {
		double i_centre = conf->mesh[i].centre;
		int j;

		k->buf[i*k->cols + i] = (1 + conf->eps_r) * 0.5;

		/* from strip[0] */
		for (j = conf->index[ID_STRIP0_START];
			j < conf->index[ID_STRIP0_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			if (j == i) {
				a->buf[i*(a->cols) + j] =
					-CONST_INV_2_PI_EPS0 * 2.0
					* conf->mesh[j].hw
					* (log(conf->mesh[j].hw) - 1);
				continue;
			}

			/* doing integral */
			sum = 0;
			y = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				if (distance < 0) {
					distance = -distance;
				}
				y = log(distance);

				sum += y;
			}

			distance = x_start - i_centre;
			if (distance < 0) {
				distance = -distance;
			}

			sum += 0.5*(log(distance) - y);

			a->buf[i*(a->cols) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}

		/* from dielectric[0] */
		for (j = conf->index[ID_DIELECTRIC0_START];
			j < conf->index[ID_DIELECTRIC0_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				if (distance < 0) {
					distance = -distance;
				}
				y = log(distance);

				sum += y;
			}

			distance = x_start - i_centre;
			if (distance < 0) {
				distance = -distance;
			}

			sum += 0.5*(log(distance) - y);

			a->buf[i*(a->cols) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}

		/* from strip[1] */
		for (j = conf->index[ID_STRIP1_START];
			j < conf->index[ID_STRIP1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double y_k, sum_k;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;
			sum_k = 0;
			y_k = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				distance = (distance*distance + h*h);

				y_k = h / distance;
				sum_k += y_k;

				y = log(sqrt(distance));
				sum += y;
			}

			distance = x_start - i_centre;
			distance = (distance*distance + h*h);

			sum_k += 0.5*((h/distance) - y_k);
			sum += 0.5*(log(sqrt(distance)) - y);

			k->buf[i*(k->cols) + j] = (1 - conf->eps_r)
							/(2*M_PI)*sum_k*dx;
			a->buf[i*(a->cols) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}

		/* from dielectric[1] */
		for (j = conf->index[ID_DIELECTRIC1_START];
			j < conf->index[ID_DIELECTRIC1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double y_k, sum_k;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;
			sum_k = 0;
			y_k = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				distance = (distance*distance + h*h);

				y_k = h / distance;
				sum_k += y_k;

				y = log(sqrt(distance));
				sum += y;
			}

			distance = x_start - i_centre;
			distance = (distance*distance + h*h);

			sum_k += 0.5*((h/distance) - y_k);
			sum += 0.5*(log(sqrt(distance)) - y);

			k->buf[i*(k->cols) + j] = (1 - conf->eps_r)
							/(2*M_PI)*sum_k*dx;
			a->buf[i*(a->cols) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}
	}
}


static INLINE void potential_strip1(
		struct MeshConfig *conf,
		struct MD *a,
		struct MD *k)

{
	double h = conf->h;
	int i;

	/* for every segments in strip[0]*/
	for (i = conf->index[ID_STRIP1_START];
		i < conf->index[ID_STRIP1_END]; ++i) {
		double i_centre = conf->mesh[i].centre;
		int j;

		k->buf[i*k->cols + i] = (1 + conf->eps_r) * 0.5;

		/* from strip[0] */
		for (j = conf->index[ID_STRIP0_START];
			j < conf->index[ID_STRIP0_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double y_k, sum_k;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;
			sum_k = 0;
			y_k = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				distance = (distance*distance + h*h);

				y_k = h / distance;
				sum_k += y_k;

				y = log(sqrt(distance));
				sum += y;
			}

			distance = x_start - i_centre;
			distance = sqrt(distance*distance + h*h);

			sum_k += 0.5*((h/distance) - y_k);
			sum += 0.5*(log(distance) - y);

			k->buf[i*(k->cols) + j] = (1 - conf->eps_r)
							/(2*M_PI)*sum_k*dx;
			a->buf[i*(a->cols) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}

		/* from dielectric[0] */
		for (j = conf->index[ID_DIELECTRIC0_START];
			j < conf->index[ID_DIELECTRIC0_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double y_k, sum_k;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;
			sum_k = 0;
			y_k = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				distance = (distance*distance + h*h);

				y_k = h / distance;
				sum_k += y_k;

				y = log(sqrt(distance));
				sum += y;
			}

			distance = x_start - i_centre;
			distance = sqrt(distance*distance + h*h);

			sum_k += 0.5*((h/distance) - y_k);
			sum += 0.5*(log(distance) - y);

			k->buf[i*(k->cols) + j] = (1 - conf->eps_r)
							/(2*M_PI)*sum_k*dx;
			a->buf[i*(a->cols) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}

		/* from strip[1] */
		for (j = conf->index[ID_STRIP1_START];
			j < conf->index[ID_STRIP1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			if (j == i) {
				a->buf[i*(a->cols) + j] =
					-CONST_INV_2_PI_EPS0 * 2.0
					* conf->mesh[j].hw
					* (log(conf->mesh[j].hw) - 1);
				continue;
			}

			/* doing integral */
			sum = 0;
			y = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				if (distance < 0) {
					distance = -distance;
				}
				y = log(distance);

				sum += y;
			}

			distance = x_start - i_centre;
			if (distance < 0) {
				distance = -distance;
			}

			sum += 0.5*(log(distance) - y);

			a->buf[i*(a->cols) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}

		/* from dielectric[1] */
		for (j = conf->index[ID_DIELECTRIC1_START];
			j < conf->index[ID_DIELECTRIC1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				if (distance < 0) {
					distance = -distance;
				}
				y = log(distance);

				sum += y;
			}

			distance = x_start - i_centre;
			if (distance < 0) {
				distance = -distance;
			}

			sum += 0.5*(log(distance) - y);

			a->buf[i*(a->cols) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}
	}
}


static INLINE void potential_dielectric0(
		struct MeshConfig *conf,
		struct MD *a,
		struct MD *k)
{
	double h = conf->h;
	int i;

	/* for every segments in dielectric[0]*/
	for (i = conf->index[ID_DIELECTRIC0_START];
		i < conf->index[ID_DIELECTRIC0_END]; ++i) {
		double i_centre = conf->mesh[i].centre;
		int j;

		k->buf[i*k->cols + i] = (1 + conf->eps_r) * 0.5;

		/* self */
		a->buf[i*(a->cols) + i] = M_PI*(1+conf->eps_r)*(1-conf->eps_r);

		/* from strip[1] */
		for (j = conf->index[ID_STRIP1_START];
			j < conf->index[ID_STRIP1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				y = h / (distance*distance + h*h);

				sum += y;
			}

			distance = x_start - i_centre;
			distance = h / (distance*distance + h*h);
			sum += 0.5*(distance - y);

			k->buf[i*(k->cols) + j] = (1 - conf->eps_r)
							/(2*M_PI)*sum*dx;
			a->buf[i*(a->cols) + j] = sum*dx;
		}

		/* from dielectric[1] */
		for (j = conf->index[ID_DIELECTRIC1_START];
			j < conf->index[ID_DIELECTRIC1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				y = h / (distance*distance + h*h);

				sum += y;
			}

			distance = x_start - i_centre;
			distance = h / (distance*distance + h*h);

			sum += 0.5*(distance - y);

			k->buf[i*(k->cols) + j] = (1 - conf->eps_r)
							/(2*M_PI)*sum*dx;
			a->buf[i*(a->cols) + j] = sum*dx;
		}
	}
}


static INLINE void potential_dielectric1(
		struct MeshConfig *conf,
		struct MD *a,
		struct MD *k)
{
	double h = conf->h;
	int i;

	/* for every segments in dielectric[1]*/
	for (i = conf->index[ID_DIELECTRIC1_START];
		i < conf->index[ID_DIELECTRIC1_END]; ++i) {
		double i_centre = conf->mesh[i].centre;
		int j;

		k->buf[i*k->cols + i] = (1 + conf->eps_r) * 0.5;

		/* self */
		a->buf[i*(a->cols) + i] = -M_PI*(1+conf->eps_r)*(1-conf->eps_r);

		/* from strip[0] */
		for (j = conf->index[ID_STRIP0_START];
			j < conf->index[ID_STRIP0_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				y = h / (distance*distance + h*h);

				sum += y;
			}

			distance = x_start - i_centre;
			distance = h / (distance*distance + h*h);
			sum += 0.5*(distance - y);

			k->buf[i*(k->cols) + j] = (1 - conf->eps_r)
							/(2*M_PI)*(sum*dx);
			a->buf[i*(a->cols) + j] = -sum*dx;
		}

		/* from dielectric[0] */
		for (j = conf->index[ID_DIELECTRIC0_START];
			j < conf->index[ID_DIELECTRIC0_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			/* doing integral */
			sum = 0;
			y = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x_start = conf->mesh[j].centre - conf->mesh[j].hw;

			for (x = x_start + dx; x <= x_end; x += dx) {
				distance = x - i_centre;
				y = h / (distance*distance + h*h);

				sum += y;
			}

			distance = x_start - i_centre;
			distance = h / (distance*distance + h*h);

			sum += 0.5*(distance - y);

			k->buf[i*(k->cols) + j] = (1 - conf->eps_r)
							/(2*M_PI)*(sum*dx);
			a->buf[i*(a->cols) + j] = -sum*dx;
		}
	}
}


static INLINE void potential_equations(
		struct MeshConfig *conf,
		struct MD *a,
		struct MD *k)
{
	potential_strip0(conf, a, k);
	potential_strip1(conf, a, k);
	potential_dielectric0(conf, a, k);
	potential_dielectric1(conf, a, k);
}


/* the factor matrix k =  (free Charge) / (all charge) at each segment */
static INLINE struct MD* mom_matrix_new(
		struct MeshConfig *conf,
		struct MD** p_k)
{
	struct MD *a, *k;

	a = md_init(conf->index[ID_MESH_CELLS],
		conf->index[ID_MESH_CELLS]);
	/* set 0 for the E equations,
	 * they are lightly sparse */
	md_fill(a, 0);

	k = md_init(conf->index[ID_MESH_CELLS],
		conf->index[ID_MESH_CELLS]);
	md_fill(k, 0);

	potential_equations(conf, a, k);

	*p_k = k;
	return a;
}

static INLINE void calc_pot(
	struct MeshConfig *conf,
	struct MD *inv_a,
	struct MD *k,
	double pot[/*2*/]
	)
{
	int i;
	struct MD *aux;
	double k1 = 0, k2 = 0;

	aux = md_mul(k, inv_a);

	for (i = conf->index[ID_STRIP_START];
		i < conf->index[ID_STRIP_END]; ++i) {
		int j;
		double hw = conf->mesh[i].hw;
		for (j = conf->index[ID_STRIP0_START];
			j < conf->index[ID_STRIP0_END]; ++j) {
			k1 -= aux->buf[i*aux->cols + j] * hw;
		}
		for (j = conf->index[ID_STRIP1_START];
			j < conf->index[ID_STRIP1_END]; ++j) {
			k2 += aux->buf[i*aux->cols + j] * hw;
		}
	}

	pot[0] = k2;
	pot[1] = k1;

	/*
	if (fabs(k1) < 1e-25) {
		pot[0] = 1;
		pot[1] = 0;
	} else {
		pot[0] = 1;
		pot[1] = k2/k1;
	}
	*/

	md_free(aux);
}


static INLINE struct MD* calc_b(
	struct MeshConfig *conf,
	struct MD *inv_a,
	double pot[/*2*/])
{
	int i;
	struct MD *b;

	b = md_init(inv_a->cols, 1);
	md_fill(b, 0);

	for (i = conf->index[ID_STRIP0_START];
		i < conf->index[ID_STRIP0_END]; ++i) {
		b->buf[i] = pot[0];
	}

	for (i = conf->index[ID_STRIP1_START];
		i < conf->index[ID_STRIP1_END]; ++i) {
		b->buf[i] = pot[1];
	}

	return b;
}

static INLINE struct MD* extract_freespace(
	struct MeshConfig *conf,
	struct MD *m)
{
	int i;
	struct MD *m0;

	m0 = md_init(conf->index[ID_STRIP_END], conf->index[ID_STRIP_END]);
	for (i = 0; i < m0->rows; ++i) {
		int j;
		for (j = 0; j < m0->cols; ++j) {
			m0->buf[i*m0->cols + j] = m->buf[i*m->cols + j];
		}
	}

	return m0;
}

static INLINE void calc_charge(
	struct MeshConfig *conf,
	struct MD *x,
	double q[]
	)
{
	int i;
	q[0] = 0;
	q[1] = 0;

	for (i = conf->index[ID_STRIP0_START];
		i < conf->index[ID_STRIP0_END]; ++i) {
		q[0] += conf->mesh[i].hw * x->buf[i];
	}

	for (i = conf->index[ID_STRIP1_START];
		i < conf->index[ID_STRIP1_END]; ++i) {
		q[1] += conf->mesh[i].hw * x->buf[i];
	}
	q[0] *= 2;
	q[1] *= 2;
}

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
	)
{
	double c[2];
	double q[2][2];
	struct MD *a[2], *k[2], *b[2], *x[2], *x_free[2];
	double pot[2][2];

	/*
	   all the arrays except `cdfs` & `cdd`
	   with the first index == 0, is in freespace
	   with the first index == 1, has a dielectric
	*/

	if (mom_trace) {
		mom_trace("Filling matrix");
	}

	a[1] = mom_matrix_new(conf, &k[1]);

	a[0] = extract_freespace(conf, a[1]);
	k[0] = md_eye(conf->index[ID_STRIP_END]);

	if (mom_trace) {
		mom_trace("Inverting matrix");
	}

	md_inverse_direct(a[0]);  /* a is inversed */
	md_inverse_direct(a[1]);

	if (mom_trace) {
		mom_trace("Extracting the result");
	}

	calc_pot(conf, a[0], k[0], pot[0]);
	calc_pot(conf, a[1], k[1], pot[1]);

	b[0] = calc_b(conf, a[0], pot[0]);
	b[1] = calc_b(conf, a[1], pot[1]);

	x[0] = md_mul(a[0], b[0]);
	x[1] = md_mul(a[1], b[1]);

	x_free[0] = md_mul(k[0], x[0]);
	x_free[1] = md_mul(k[1], x[1]);

	calc_charge(conf, x_free[0], q[0]);
	calc_charge(conf, x_free[1], q[1]);

	md_free(a[0]);
	md_free(a[1]);
	md_free(k[0]);
	md_free(k[1]);
	md_free(b[0]);
	md_free(b[1]);

	if (cdfs) {
		cdfs[0] = x_free[0];
		cdfs[1] = x[0];
	} else {
		md_free(x[0]);
		md_free(x_free[0]);
	}

	if (cdd) {
		cdd[0] = x_free[1];
		cdd[1] = x[1];
	} else {
		md_free(x[1]);
		md_free(x_free[1]);
	}

	c[0] = q[0][0] / (pot[0][0] - pot[0][1]);
	c[1] = q[1][0] / (pot[1][0] - pot[1][1]);

	if (capacity) {
		capacity[0] = c[0];
		capacity[1] = c[1];
	}

	return 1/(C0*sqrt(c[0]*c[1]));
}
