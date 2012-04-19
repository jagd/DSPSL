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
	int h = conf->h;
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
			double x, dx, x_end, sum, y_old;
			double distance;

			if (j == i) {
				a->buf[i*(a->rows) + j] = 
					CONST_INV_2_PI_EPS0 * 2.0
					* conf->mesh[j].hw
					* (log(conf->mesh[j].hw) - 1);
				continue;
			}

			/* doing integral */
			sum = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x = conf->mesh[j].centre - conf->mesh[j].hw;

			distance = x - i_centre;
			if (distance <= 0) {
				distance = -distance;
			}

			y_old = log(distance);
			for (x = x + dx; x < x_end; x += dx) {
				int y;

				distance = x - i_centre;
				if (distance <= 0) {
					distance = -distance;
				} y = log(distance);

				sum += (y + y_old) * dx;
				y_old = y;
			}

			a->buf[i*(a->rows) + j] = sum * -0.5;
		}

		/* from dielectric[0] */
		for (j = conf->index[ID_DIELECTRIC0_START];
			j < conf->index[ID_DIELECTRIC0_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_end, sum, y_old;
			double distance;

			/* doing integral */
			sum = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x = conf->mesh[j].centre - conf->mesh[j].hw;

			distance = x - i_centre;
			if (distance <= 0) {
				distance = -distance;
			}

			y_old = log(distance);
			for (x = x + dx; x < x_end; x += dx) {
				int y;

				distance = x - i_centre;
				if (distance <= 0) {
					distance = -distance;
				}
				y = log(distance);

				sum += (y + y_old) * dx;
				y_old = y;
			}

			a->buf[i*(a->rows) + j] = sum * -0.5;
		}

		/* from strip[1] */
		for (j = conf->index[ID_STRIP1_START];
			j < conf->index[ID_STRIP1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_end, sum, y_old;
			double distance;

			/* doing integral */
			sum = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x = conf->mesh[j].centre - conf->mesh[j].hw;

			distance = x - i_centre;
			distance = sqrt(distance*distance + h*h);

			y_old = log(distance);
			for (x = x + dx; x < x_end; x += dx) {
				int y;

				distance = x - i_centre;
				distance = sqrt(distance*distance + h*h);
				y = log(distance);

				sum += (y + y_old) * dx;
				y_old = y;
			}

			a->buf[i*(a->rows) + j] = sum * -0.5;
		}

		/* from dielectric[1] */
		for (j = conf->index[ID_DIELECTRIC1_START];
			j < conf->index[ID_DIELECTRIC1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_end, sum, y_old;
			double distance;

			/* doing integral */
			sum = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x = conf->mesh[j].centre - conf->mesh[j].hw;

			distance = x - i_centre;
			distance = sqrt(distance*distance + h*h);

			y_old = log(distance);
			for (x = x + dx; x < x_end; x += dx) {
				int y;

				distance = x - i_centre;
				distance = sqrt(distance*distance + h*h);
				y = log(distance);

				sum += (y + y_old) * dx;
				y_old = y;
			}

			a->buf[i*(a->rows) + j] = sum * -0.5;
		}
	}
}

static INLINE void potential_strip1(struct MeshConfig *conf, struct MD *a)
{
	int h = conf->h;
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
			double x, dx, x_end, sum, y_old;
			double distance;

			/* doing integral */
			sum = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x = conf->mesh[j].centre - conf->mesh[j].hw;

			distance = x - i_centre;
			distance = sqrt(distance*distance + h*h);

			y_old = log(distance);
			for (x = x + dx; x < x_end; x += dx) {
				int y;

				distance = x - i_centre;
				distance = sqrt(distance*distance + h*h);
				y = log(distance);

				sum += (y + y_old) * dx;
				y_old = y;
			}

			a->buf[i*(a->rows) + j] = sum * -0.5;
		}

		/* from dielectric[0] */
		for (j = conf->index[ID_DIELECTRIC0_START];
			j < conf->index[ID_DIELECTRIC0_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_end, sum, y_old;
			double distance;

			/* doing integral */
			sum = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x = conf->mesh[j].centre - conf->mesh[j].hw;

			distance = x - i_centre;
			distance = sqrt(distance*distance + h*h);

			y_old = log(distance);
			for (x = x + dx; x < x_end; x += dx) {
				int y;

				distance = x - i_centre;
				distance = sqrt(distance*distance + h*h);
				y = log(distance);

				sum += (y + y_old) * dx;
				y_old = y;
			}

			a->buf[i*(a->rows) + j] = sum * -0.5;
		}

		/* from strip[1] */
		for (j = conf->index[ID_STRIP1_START];
			j < conf->index[ID_STRIP1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_end, sum, y_old;
			double distance;

			if (j == i) {
				a->buf[i*(a->rows) + j] = 
					CONST_INV_2_PI_EPS0 * 2.0
					* conf->mesh[j].hw
					* (log(conf->mesh[j].hw) - 1);
				continue;
			}

			/* doing integral */
			sum = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x = conf->mesh[j].centre - conf->mesh[j].hw;

			distance = x - i_centre;
			if (distance <= 0) {
				distance = -distance;
			}

			y_old = log(distance);
			for (x = x + dx; x < x_end; x += dx) {
				int y;

				distance = x - i_centre;
				if (distance <= 0) {
					distance = -distance;
				} y = log(distance);

				sum += (y + y_old) * dx;
				y_old = y;
			}

			a->buf[i*(a->rows) + j] = sum * -0.5;
		}

		/* from dielectric[1] */
		for (j = conf->index[ID_DIELECTRIC1_START];
			j < conf->index[ID_DIELECTRIC1_END]; ++j) {

			/* variables for the integral */
			/* `dx` always > 0*/
			double x, dx, x_end, sum, y_old;
			double distance;

			/* doing integral */
			sum = 0;

			dx = conf->mesh[j].hw * 2.0 / INTEGRAL_STEPS;
			x_end = conf->mesh[j].centre + conf->mesh[j].hw;
			x = conf->mesh[j].centre - conf->mesh[j].hw;

			distance = x - i_centre;
			if (distance <= 0) {
				distance = -distance;
			}

			y_old = log(distance);
			for (x = x + dx; x < x_end; x += dx) {
				int y;

				distance = x - i_centre;
				if (distance <= 0) {
					distance = -distance;
				}
				y = log(distance);

				sum += (y + y_old) * dx;
				y_old = y;
			}

			a->buf[i*(a->rows) + j] = sum * -0.5;
		}
	}
}

static INLINE void potential_equations(struct MeshConfig *conf, struct MD *a)
{
	potential_strip0(conf, a);
	potential_strip1(conf, a);
}

struct MD* mom_matrix_new(struct MeshConfig *conf)
{
	struct MD *a;

	a = md_init(conf->index[ID_MESH_CELLS],
		conf->index[ID_MESH_CELLS]);

md_fill(a, 0);

	potential_equations(conf, a);

md_dump(a);
	// TODO: 
	return a;
}
