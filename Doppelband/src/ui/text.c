#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../mom_mesh.h"
#include "../mom.h"

#define PATH_GNUPLOT "/usr/bin/gnuplot"

#define PATH_LENGTH 2000

double w1, w2, d, h, eps_r, port_ext;
char path[PATH_LENGTH];

#define DEFAULT_PLOT_PATH "plot.txt"

static void chomp(char str[])
{
	int i;

	i = 0;
	while (str[i]) {
		++i;
	}
	--i;

	while (i >= 0 &&
		(str[i] == '\n'
		 || str[i] == '\r'
		 || str[i] == ' '
		 || str[i] == '\t'))
	{
		str[i--] = 0;
	}
}

static void input()
{
	puts(
		"\n\n"
		"                              w1   |< extra width >|\n"
		"    ------------------------=======-----------------\n"
		"    ////////////////////////////////////////////////\n"
		"    ////////////////////////////////////////////////\n"
		"    /////////////////    epsilon   /////////////////\n"
		"    ////////////////////////////////////////////////\n"
		"    /////////////////   w2  ////////////////////////\n"
		"    -----------------=========----------------------\n"
		"    |< extra width >|    |      |\n"
		"                         |< d  >|\n"
		"\n\n"
	      );

	while (1) {
		printf("Width `w1` of the top strip (in mm) = ");
		scanf("%le", &w1);
		w1 *= 1e-3;
		if (w1 < 1e-20) {
			puts("`w1` must be bigger than zero!");
		} else {
			break;
		}
	}

	while (1) {
		printf("Width `w2` of the bottom strip (in mm) = ");
		scanf("%le", &w2);
		w2 *= 1e-3;

		if (w2 < 1e-20) {
			puts("`w2` must be bigger than zero!");
		} else {
			break;
		}
	}

	while (1) {
		printf("Extra width of the port (in mm) = ");
		scanf("%le", &port_ext);
		port_ext *= 1e-3;

		if (port_ext < 1e-20) {
			puts("The extra width must be bigger than zero!");
		} else {
			break;
		}
	}

	while (1) {
		printf("Height of the substrat (in mm) = ");
		scanf("%le", &h);
		h *= 1e-3;

		if (h < 1e-20) {
			puts("The Height must be bigger than zero!");
		} else {
			break;
		}
	}

	printf("The offset `d` between the centre of both strips (in mm) = ");
	scanf("%le", &d);
	d *= 1e-3;

	while (1) {
		printf("Relative permittivity = ");
		scanf("%le", &eps_r);

		if (eps_r < 1.2) {
			puts("Too small permittivity is not allowed!");
		} else {
			break;
		}
	}

	fgets(path, PATH_LENGTH - 1, stdin); /* pass one \n */

	printf("The file name to plot the charge density with dielectric "
		"[%s] : ", DEFAULT_PLOT_PATH);
	fgets(path, PATH_LENGTH - 1, stdin);
	chomp(path);
	if (strlen(path) == 0) {
		strcpy(path, DEFAULT_PLOT_PATH);
	}
}

void calc()
{
	struct MeshConfig *conf;
	struct MD *x[2], *x0[2];
	int i;
	double c[2];
	double z0;

	FILE *fplot;

	conf = mesh_new(
		w1, w2,
		d,
		port_ext,
		h,
		eps_r
		);

	printf("Mesh cells = %d\n", conf->index[ID_MESH_CELLS]);
	printf("\tcells for strips = %d\n", conf->index[ID_STRIP_END]);
	printf("\tcells for dielectrics = %d\n",
			conf->index[ID_DIELECTRIC_END]
				-conf->index[ID_DIELECTRIC_START]);

	z0 = mom(conf, x0, x, c);


	printf("Effective permittivity = %lf\n", c[1]/c[0]);
	printf("Z0 = %lf Ohm\n", z0);

	fplot = fopen(path, "w");

	if (fplot) {
		fprintf(fplot, "#!"PATH_GNUPLOT"\n");

		fprintf(fplot, "set title \"all charges with dielectric\"\n");
		fprintf(fplot, "plot '-' notitle with impulse, \\\n"
				"\t'-' notitle with impulse\n");
		fprintf(fplot, "# top all charges\n");
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


		fprintf(fplot, "set title \"free charges with dielectric\"\n");
		fprintf(fplot, "plot '-' notitle with impulse, \\\n"
				"\t'-' notitle with impulse\n");
		fprintf(fplot, "# top free charges\n");
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


		fprintf(fplot, "set title \"charges in free space\"\n");
		fprintf(fplot, "plot '-' notitle with impulse, \\\n"
				"\t'-' notitle with impulse\n");
		fprintf(fplot, "# top free charges\n");
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

	input();
	calc();

	return 0;
}
