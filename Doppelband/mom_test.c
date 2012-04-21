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
	/*
	struct MD *x, *x_free;
	int i;
	*/
	double c[2];
	double z0;

	conf = mesh_new(
		2.4e-3, 10e-3,
		0,
		3e-3,
		0.79e-3,
		2.2
		);

	z0 = mom(conf, NULL, NULL, c);

	printf("C0 = %le F        C1 = %le F\n", c[0], c[1]);
	printf("Z0 = %lf Ohm\n", z0);

	/*

	printf("# top all charges\n");

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

	printf("# bottom all charges\n");

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

	printf("# top free charges\n");

	for (i = conf->index[ID_STRIP0_START];
		i < conf->index[ID_STRIP0_END]; ++i) {
		printf("%le %le %le\n", conf->mesh[i].centre, x_free->buf[i]
				, conf->mesh[i].hw);
	}

	for (i = conf->index[ID_DIELECTRIC0_START];
		i < conf->index[ID_DIELECTRIC0_END]; ++i) {
		printf("%le %le %le\n", conf->mesh[i].centre, x_free->buf[i]
				, conf->mesh[i].hw);
	}

	printf("# bottom free charges\n");

	for (i = conf->index[ID_STRIP1_START];
		i < conf->index[ID_STRIP1_END]; ++i) {
		printf("%le %le %le\n", conf->mesh[i].centre, x_free->buf[i]
				, conf->mesh[i].hw);
	}

	for (i = conf->index[ID_DIELECTRIC1_START];
		i < conf->index[ID_DIELECTRIC1_END]; ++i) {
		printf("%le %le %le\n", conf->mesh[i].centre, x_free->buf[i]
				, conf->mesh[i].hw);
	}



	md_free(x);
	md_free(x_free);
	*/
	mesh_free(conf);
}

int main()
{
	test_matrix();
	return 0;
}
