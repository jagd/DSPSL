#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../global.h"
#include "../mom_mesh.h"
#include "../mom.h"

#define PATH_GNUPLOT TEXT("/usr/bin/gnuplot")

#define PATH_LENGTH 2000

double w1, w2, d, h, eps_r, port_ext;
char path[PATH_LENGTH];

#define DEFAULT_PLOT_PATH TEXT("plot.txt")

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
		TEXT("\n\n")
		TEXT("                              w1   |< extra width >|\n")
		TEXT("    ------------------------=======-----------------\n")
		TEXT("    ////////////////////////////////////////////////\n")
		TEXT("    ////////////////////////////////////////////////\n")
		TEXT("    /////////////////    epsilon   /////////////////\n")
		TEXT("    ////////////////////////////////////////////////\n")
		TEXT("    /////////////////   w2  ////////////////////////\n")
		TEXT("    -----------------=========----------------------\n")
		TEXT("    |< extra width >|    |      |\n")
		TEXT("                         |< d  >|\n")
		TEXT("\n\n")
	      );

	while (1) {
		printf(TEXT("Width `w1` of the top strip (in mm) = "));
		scanf(TEXT("%le"), &w1);
		w1 *= 1e-3;
		if (w1 < 1e-20) {
			puts(TEXT("`w1` must be bigger than zero!"));
		} else {
			break;
		}
	}

	while (1) {
		printf(TEXT("Width `w2` of the bottom strip (in mm) = "));
		scanf(TEXT("%le"), &w2);
		w2 *= 1e-3;

		if (w2 < 1e-20) {
			puts(TEXT("`w2` must be bigger than zero!"));
		} else {
			break;
		}
	}

	while (1) {
		printf(TEXT("Extra width of the port (in mm) = "));
		scanf(TEXT("%le"), &port_ext);
		port_ext *= 1e-3;

		if (port_ext < 1e-20) {
			puts(TEXT("The extra width must be bigger than zero!"));
		} else {
			break;
		}
	}

	while (1) {
		printf(TEXT("Height of the substrat (in mm) = "));
		scanf(TEXT("%le"), &h);
		h *= 1e-3;

		if (h < 1e-20) {
			puts(TEXT("The Height must be bigger than zero!"));
		} else {
			break;
		}
	}

	printf(TEXT("The offset `d` between the centre of both strips (in mm) = "));
	scanf(TEXT("%le"), &d);
	d *= 1e-3;

	while (1) {
		printf(TEXT("Relative permittivity = "));
		scanf(TEXT("%le"), &eps_r);

		if (eps_r < 1) {
			puts(TEXT("The relative permittivity must > 1.0 !"));
		} else {
			break;
		}
	}

	fgets(path, PATH_LENGTH - 1, stdin); /* pass one \n */

	printf(TEXT("The file name to plot the charge density with dielectric ")
		TEXT("[%s] : "), DEFAULT_PLOT_PATH);
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

	printf(TEXT("Mesh cells = %d\n"), conf->index[ID_MESH_CELLS]);
	printf(TEXT("\tcells for strips = %d\n"), conf->index[ID_STRIP_END]);
	printf(TEXT("\tcells for dielectrics = %d\n"),
			conf->index[ID_DIELECTRIC_END]
				-conf->index[ID_DIELECTRIC_START]);

	z0 = mom(conf, x0, x, c);


	printf(TEXT("Effective permittivity = %lf\n"), c[1]/c[0]);
	printf(TEXT("Z0 = %lf Ohm\n"), z0);

	fplot = fopen(path, TEXT("w"));

	if (fplot) {
		fprintf(fplot, TEXT("#!") PATH_GNUPLOT TEXT("\n"));

		fprintf(fplot, TEXT("set title \"all charges with dielectric\"\n"));
		fprintf(fplot, TEXT("plot '-' notitle with impulse, \\\n")
				TEXT("\t'-' notitle with impulse\n"));
		fprintf(fplot, TEXT("# top all charges\n"));
		for (i = conf->index[ID_STRIP0_START];
			i < conf->index[ID_STRIP0_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[1]->buf[i]);
		}
		for (i = conf->index[ID_DIELECTRIC0_START];
			i < conf->index[ID_DIELECTRIC0_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[1]->buf[i]);
		}
		fprintf(fplot, TEXT("e\n"));
		fprintf(fplot, TEXT("# bottom all charges\n"));
		for (i = conf->index[ID_STRIP1_START];
			i < conf->index[ID_STRIP1_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[1]->buf[i]);
		}
		for (i = conf->index[ID_DIELECTRIC1_START];
			i < conf->index[ID_DIELECTRIC1_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[1]->buf[i]);
		}
		fprintf(fplot, TEXT("e\npause -1\n"));


		fprintf(fplot, TEXT("set title \"free charges with dielectric\"\n"));
		fprintf(fplot, TEXT("plot '-' notitle with impulse, \\\n")
				TEXT("\t'-' notitle with impulse\n"));
		fprintf(fplot, TEXT("# top free charges\n"));
		for (i = conf->index[ID_STRIP0_START];
			i < conf->index[ID_STRIP0_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[0]->buf[i]);
		}
		for (i = conf->index[ID_DIELECTRIC0_START];
			i < conf->index[ID_DIELECTRIC0_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[0]->buf[i]);
		}
		fprintf(fplot, TEXT("e\n"));
		fprintf(fplot, TEXT("# bottom free charges\n"));
		for (i = conf->index[ID_STRIP1_START];
			i < conf->index[ID_STRIP1_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[0]->buf[i]);
		}
		for (i = conf->index[ID_DIELECTRIC1_START];
			i < conf->index[ID_DIELECTRIC1_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[0]->buf[i]);
		}
		fprintf(fplot, TEXT("e\npause -1\n"));


		fprintf(fplot, TEXT("set title \"charges in free space\"\n"));
		fprintf(fplot, TEXT("plot '-' notitle with impulse, \\\n")
				TEXT("\t'-' notitle with impulse\n"));
		fprintf(fplot, TEXT("# top free charges\n"));
		for (i = conf->index[ID_STRIP0_START];
			i < conf->index[ID_STRIP0_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[1]->buf[i]);
		}
		fprintf(fplot, TEXT("e\n"));
		fprintf(fplot, TEXT("# bottom free charges\n"));
		for (i = conf->index[ID_STRIP1_START];
			i < conf->index[ID_STRIP1_END]; ++i) {
			fprintf(fplot, TEXT("%le %le\n"), conf->mesh[i].centre, x[1]->buf[i]);
		}
		fprintf(fplot, TEXT("e\npause -1\n"));

		fclose(fplot);
	}

	md_free(x[0]);
	md_free(x[1]);
	md_free(x0[0]);
	md_free(x0[1]);

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
	mom_error=mom_trace;

	input();
	calc();

	return 0;
}
