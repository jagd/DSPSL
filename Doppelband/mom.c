#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mom_mesh.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define C0	3e8
#define Z0	(120 * M_PI)
#define MU0	(Z0 / C0)
#define EPS0	(1 / (Z0*C0))


void TestCase()
{
	struct MeshConfig *conf;
	int i;

	conf = mesh_new(
		2.4e-3, 5e-3,
		0,
		3e-3,
		0.79e-3
		);

	puts("#Strip[0]:");
	for (i = conf->index[ID_STRIP0_START];
		i < conf->index[ID_STRIP0_END]; ++i) {
		printf("%lf %lf 2\n", conf->mesh[i].centre, conf->mesh[i].hw);
	}

	puts("#Dielectric[0]:");
	for (i = conf->index[ID_DIELECTRIC0_START];
		i < conf->index[ID_DIELECTRIC0_END]; ++i) {
		printf("%lf %lf 1\n", conf->mesh[i].centre, conf->mesh[i].hw);
	}


	puts("#Strip[1]:");
	for (i = conf->index[ID_STRIP1_START];
		i < conf->index[ID_STRIP1_END]; ++i) {
		printf("%lf %lf 2\n", conf->mesh[i].centre, conf->mesh[i].hw);
	}

	puts("#Dielectric[1]:");
	for (i = conf->index[ID_DIELECTRIC1_START];
		i < conf->index[ID_DIELECTRIC1_END]; ++i) {
		printf("%lf %lf 1\n", conf->mesh[i].centre, conf->mesh[i].hw);
	}

	mesh_free(conf);
}


int main()
{
	TestCase();

	return 0;
}
