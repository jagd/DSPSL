#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mom_mesh.h"
#include "mom.h"

#define PATH_GNUPLOT "/usr/bin/gnuplot"

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

	mesh_free(conf);
}

void test_matrix()
{
	struct MeshConfig *conf;
	struct MD *x[2], *x0[2];
	int i;
	double c[2];
	double z0;

	FILE *fplot;

	conf = mesh_new(
		2.4e-3, 10e-3,
		0,
		3e-3,
		0.79e-3,
		2.2
		);

	z0 = mom(conf, x0, x, c);


	printf("C0 = %le F        C1 = %le F\n", c[0], c[1]);
	printf("Z0 = %lf Ohm\n", z0);

	fplot = fopen("dielectric.plot", "w");

	if (fplot) {
		fprintf(fplot, "#!"PATH_GNUPLOT"\n");
		fprintf(fplot, "# top all charges\n");
		fprintf(fplot, "plot '-' notitle with impulse, \\\n"
				"\t'-' notitle with impulse\n");

		for (i = conf->index[ID_STRIP0_START];
			i < conf->index[ID_STRIP0_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[1]->buf[i]);
		}

		for (i = conf->index[ID_DIELECTRIC0_START];
			i < conf->index[ID_DIELECTRIC0_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[1]->buf[i]);
		}
		fprintf(fplot, "e\n");

		fprintf(fplot, "# bottom all charges\n");

		for (i = conf->index[ID_STRIP1_START];
			i < conf->index[ID_STRIP1_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[1]->buf[i]);
		}

		for (i = conf->index[ID_DIELECTRIC1_START];
			i < conf->index[ID_DIELECTRIC1_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[1]->buf[i]);
		}
		fprintf(fplot, "e\npause -1\n");
		/*****/
		fprintf(fplot, "# top free charges\n");
		fprintf(fplot, "plot '-' notitle with impulse, \\\n"
				"\t'-' notitle with impulse\n");

		for (i = conf->index[ID_STRIP0_START];
			i < conf->index[ID_STRIP0_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[0]->buf[i]);
		}

		for (i = conf->index[ID_DIELECTRIC0_START];
			i < conf->index[ID_DIELECTRIC0_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[0]->buf[i]);
		}
		fprintf(fplot, "e\n");

		fprintf(fplot, "# bottom free charges\n");

		for (i = conf->index[ID_STRIP1_START];
			i < conf->index[ID_STRIP1_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[0]->buf[i]);
		}

		for (i = conf->index[ID_DIELECTRIC1_START];
			i < conf->index[ID_DIELECTRIC1_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[0]->buf[i]);
		}
		fprintf(fplot, "e\npause -1\n");
		fclose(fplot);
	}

	fplot = fopen("freespace.plot", "w");

	if (fplot) {
		fprintf(fplot, "#!"PATH_GNUPLOT"\n");
		fprintf(fplot, "# top free charges\n");
		fprintf(fplot, "plot '-' notitle with impulse, \\\n"
				"\t'-' notitle with impulse\n");

		for (i = conf->index[ID_STRIP0_START];
			i < conf->index[ID_STRIP0_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[1]->buf[i]);
		}

		fprintf(fplot, "e\n");

		fprintf(fplot, "# bottom free charges\n");

		for (i = conf->index[ID_STRIP1_START];
			i < conf->index[ID_STRIP1_END]; ++i) {
			fprintf(fplot, "%le %le\n", conf->mesh[i].centre, x[1]->buf[i]);
		}

		fprintf(fplot, "e\npause -1\n");
		fclose(fplot);
	}
	md_free(x[0]);
	md_free(x[1]);

	mesh_free(conf);
}

/*
   wrap the puts function,
   because its `const char*` parameter is not C89 compatible
*/
void trace(char *s)
{
	puts(s);
}

int main()
{
	mom_trace=trace;

	test_matrix();

	return 0;
}
