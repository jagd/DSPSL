#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mom_mesh.h"
#include "md.h"

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

static INLINE void potential_strip0(struct MeshConfig *conf, struct MD *a)
{
	double h = conf->h;
	int i;

	/* for every segments in strip[0]*/
	for (i = conf->index[ID_STRIP0_START];
		i < conf->index[ID_STRIP0_END]; ++i) {
		double i_centre = conf->mesh[i].centre;
		int j;

		/* from strip[0] */
		for (j = conf->index[ID_STRIP0_START];
			j < conf->index[ID_STRIP0_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			if (j == i) {
				a->buf[i*(a->rows) + j] =
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

			a->buf[i*(a->rows) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
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

			a->buf[i*(a->rows) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}

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
				distance = sqrt(distance*distance + h*h);
				y = log(distance);

				sum += y;
			}

			distance = x_start - i_centre;
			distance = sqrt(distance*distance + h*h);

			sum += 0.5*(log(distance) - y);

			a->buf[i*(a->rows) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
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
				distance = sqrt(distance*distance + h*h);
				y = log(distance);

				sum += y;
			}

			distance = x_start - i_centre;
			distance = sqrt(distance*distance + h*h);

			sum += 0.5*(log(distance) - y);

			a->buf[i*(a->rows) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}
	}
}


static INLINE void potential_strip1(struct MeshConfig *conf, struct MD *a)
{
	double h = conf->h;
	int i;

	/* for every segments in strip[0]*/
	for (i = conf->index[ID_STRIP1_START];
		i < conf->index[ID_STRIP1_END]; ++i) {
		double i_centre = conf->mesh[i].centre;
		int j;

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
				distance = sqrt(distance*distance + h*h);
				y = log(distance);

				sum += y;
			}

			distance = x_start - i_centre;
			distance = sqrt(distance*distance + h*h);

			sum += 0.5*(log(distance) - y);

			a->buf[i*(a->rows) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
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
				distance = sqrt(distance*distance + h*h);
				y = log(distance);

				sum += y;
			}

			distance = x_start - i_centre;
			distance = sqrt(distance*distance + h*h);

			sum += 0.5*(log(distance) - y);

			a->buf[i*(a->rows) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}

		/* from strip[1] */
		for (j = conf->index[ID_STRIP1_START];
			j < conf->index[ID_STRIP1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_start, x_end, sum, y;
			double distance;

			if (j == i) {
				a->buf[i*(a->rows) + j] =
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

			a->buf[i*(a->rows) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
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

			a->buf[i*(a->rows) + j] = -CONST_INV_2_PI_EPS0*sum*dx;
		}
	}
}


static INLINE void potential_dielectric0(struct MeshConfig *conf, struct MD *a)
{
	double h = conf->h;
	int i;

	/* for every segments in dielectric[0]*/
	for (i = conf->index[ID_DIELECTRIC0_START];
		i < conf->index[ID_DIELECTRIC0_END]; ++i) {
		double i_centre = conf->mesh[i].centre;
		int j;

		/* self */
		a->buf[i*(a->rows) + i] = M_PI*(1+conf->eps_r)*(1-conf->eps_r);

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

			a->buf[i*(a->rows) + j] = sum*dx;
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

			a->buf[i*(a->rows) + j] = sum*dx;
		}
	}
}


static INLINE void potential_dielectric1(struct MeshConfig *conf, struct MD *a)
{
	double h = conf->h;
	int i;

	/* for every segments in dielectric[1]*/
	for (i = conf->index[ID_DIELECTRIC1_START];
		i < conf->index[ID_DIELECTRIC1_END]; ++i) {
		double i_centre = conf->mesh[i].centre;
		int j;

		/* self */
		a->buf[i*(a->rows) + i] = -M_PI*(1+conf->eps_r)*(1-conf->eps_r);

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

			a->buf[i*(a->rows) + j] = -sum*dx;
		}

		/* from dielectric[1] */
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

			a->buf[i*(a->rows) + j] = -sum*dx;
		}
	}
}


static INLINE void potential_equations(struct MeshConfig *conf, struct MD *a)
{
	potential_strip0(conf, a);
	potential_strip1(conf, a);
	potential_dielectric0(conf, a);
	potential_dielectric1(conf, a);
}


struct MD* mom_matrix_new(struct MeshConfig *conf)
{
	struct MD *a;

	a = md_init(conf->index[ID_MESH_CELLS],
		conf->index[ID_MESH_CELLS]);

	/* set 0 for the E equations,
	 * they are lightly sparse */

	md_fill(a, 0);

	potential_equations(conf, a);

	return a;
}


struct MD* calc_charge(
		struct MeshConfig *conf,
		struct MD* a, /* will be dirty */
		double charge[/* 2 */] )
{
	double q[2];
	struct MD *b, *x;
	int i;
	double pot[2];

	b = md_init(a->cols, 1);
	md_fill(b, 0);

	pot[0] = 1;
	pot[1] = -1;

	for (i = conf->index[ID_STRIP0_START];
		i < conf->index[ID_STRIP0_END]; ++i) {
		b->buf[i] = pot[0];
	}

	for (i = conf->index[ID_STRIP1_START];
		i < conf->index[ID_STRIP1_END]; ++i) {
		b->buf[i] = pot[1];
	}

	md_inverse_direct(a);
	x = md_mul(a, b);

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

	charge[0] = q[0]*2;
	charge[1] = q[1]*2;

	return x;
}
