#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mom_mesh.h"
#include "mom.h"

void test_mesh()
{
	struct MeshConfig *conf;
	/* int i; */

	conf = mesh_new(
		2.4e-3, 5e-3,
		0,
		3e-3,
		0.79e-3,
		2.2
		);

	/*
	puts("#Strip[0]:");
	for (i = conf->index[ID_STRIP0_START];
		i < conf->index[ID_STRIP0_END]; ++i) {
		printf("%le %le 2\n", conf->mesh[i].centre, conf->mesh[i].hw);
	}

	puts("#Dielectric[0]:");
	for (i = conf->index[ID_DIELECTRIC0_START];
		i < conf->index[ID_DIELECTRIC0_END]; ++i) {
		printf("%le %le 1\n", conf->mesh[i].centre, conf->mesh[i].hw);
	}


	puts("#Strip[1]:");
	for (i = conf->index[ID_STRIP1_START];
		i < conf->index[ID_STRIP1_END]; ++i) {
		printf("%le %le 2\n", conf->mesh[i].centre, conf->mesh[i].hw);
	}

	puts("#Dielectric[1]:");
	for (i = conf->index[ID_DIELECTRIC1_START];
		i < conf->index[ID_DIELECTRIC1_END]; ++i) {
		printf("%le %le 1\n", conf->mesh[i].centre, conf->mesh[i].hw);
	}
	*/

	mesh_free(conf);
}

void test_matrix()
{
	struct MeshConfig *conf;
	struct MD *x;
	int i;

	double q[2];

	conf = mesh_new(
		2.4e-3, 5e-3,
		0,
		3e-3,
		0.79e-3,
		2.2
		);

	x = calc_charge(conf, q);

	printf("Q[0] = %le\tQ[1] = %le\n", q[0], q[1]);

	printf("# top\n");

	for (i = conf->index[ID_STRIP0_START];
		i < conf->index[ID_STRIP0_END]; ++i) {
		printf("%le %le %le\n", conf->mesh[i].centre, x->buf[i]
				, conf->mesh[i].hw);
	}

	for (i = conf->index[ID_DIELECTRIC0_START];
		i < conf->index[ID_DIELECTRIC0_END]; ++i) {
		printf("%le %le %le\n", conf->mesh[i].centre, x->buf[i]
				, conf->mesh[i].hw);
	}

	printf("# bottom\n");

	for (i = conf->index[ID_STRIP1_START];
		i < conf->index[ID_STRIP1_END]; ++i) {
		printf("%le %le %le\n", conf->mesh[i].centre, x->buf[i]
				, conf->mesh[i].hw);
	}

	for (i = conf->index[ID_DIELECTRIC1_START];
		i < conf->index[ID_DIELECTRIC1_END]; ++i) {
		printf("%le %le %le\n", conf->mesh[i].centre, x->buf[i]
				, conf->mesh[i].hw);
	}

	mesh_free(conf);

	md_free(x);
}

int main()
{
	test_matrix();
	return 0;
}
