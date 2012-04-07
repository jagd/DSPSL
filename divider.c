#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define Z_FS 376.73
#define C0 3e8

#define BAR() puts("\n========================================\n")
#define STR_PREFIX ">>> "

double calc_lambda(double freq, double eps_r_eff)
{
	return C0 / freq / sqrt(eps_r_eff);
}

double w2z(double h, double eps_r, double w, double *p_eps_r_eff)
{
	double buff_eps;
	double f;
	double z;

	if (p_eps_r_eff == NULL) {
		p_eps_r_eff = &buff_eps;
	}

	if (w / h < 1) {
		f = 1 / sqrt(1 + 12*h/w) + 0.04 * (1 - w/h) * (1 - w/h);
		*p_eps_r_eff = 0.5 * (eps_r + 1) + 0.5 * (eps_r - 1) * f;
		z = 60 / sqrt(*p_eps_r_eff) * log(8 * h/w + 0.25 * w/h);
	} else {
		f = 1 / sqrt(1 + 12*h/w);
		*p_eps_r_eff = 0.5 * (eps_r + 1) + 0.5 * (eps_r - 1) * f;
		z = 120 * M_PI / ( sqrt(*p_eps_r_eff)
				* (w/h + 1.393 + 0.667 * log(w/h + 1.44)) );
	}

	return z;
}

double z2w(double h, double eps_r, double z)
{
	/* Width has the same Unit as Height */

	double a, b;
	double lt2, gt2;

	a = M_PI / Z_FS * sqrt(2*(eps_r + 1)) * z
		+ (eps_r - 1) / (eps_r + 1) * (0.226 + 0.121 / eps_r);

	lt2 = 4 / (0.5 * exp(a) - exp(-a));

	if (lt2 > 0 && lt2 <= 2) {
		return lt2 * h;
	}

	/* else */

	b = M_PI * Z_FS / (2 * z * sqrt(eps_r));

	gt2 = (eps_r - 1) / (M_PI * eps_r) * (log(b - 1) + 0.293 - 0.517/eps_r)
		+ 2 / M_PI * (b - 1 - log(2*b - 1));

	if (gt2 >= 2) {
		return gt2 * h;
	}

	fprintf(stderr, "error by calculating the width\n");
	return 0;
}


void branchline(double k2, double z0, double h, double eps_r, double freq)
{
	double zh, zv;
	double wh, wv;
	double lambda_h, lambda_v;
	double tmp_eps_eff;

	printf(
		STR_PREFIX"Branchline:\n\n"
		"Port 1                     Port 2 \n"
		"  o-----+-----[Zh]-----+-----o    \n"
		"        |              |          \n"
		"        |              |          \n"
		"       [Zv]           [Zv]        \n"
		"        |              |          \n"
		"        |              |          \n"
		"    *---+-----[Zh]-----+-----o    \n"
		"    |                     Port 3 \n"
		"    ^                            \n"
		"    Z0                           \n"
		"    v                            \n"
		"    |                            \n"
		"   ---                           \n"
		"    -                            \n"
						"\n\n"
		);

	zh = z0 * sqrt(k2 / (1 + k2));
	zv = z0 * sqrt(k2);


	wh = z2w(h, eps_r, zh);
	wv = z2w(h, eps_r, zv);

	w2z(h, eps_r, wh, &tmp_eps_eff);
	lambda_h = calc_lambda(freq, tmp_eps_eff) / 1e6;
	w2z(h, eps_r, wv, &tmp_eps_eff);
	lambda_v = calc_lambda(freq, tmp_eps_eff) / 1e6;

	printf(
		STR_PREFIX"Zh = %lf Ohm ,\tW = %lf mm ,\t Length = %lf mm\n"
		STR_PREFIX"Zv = %lf Ohm ,\tW = %lf mm ,\t Length = %lf mm\n"
		, zh, wh, lambda_h/4, zv, wv, lambda_v/4);
}

void wilkinson(double k2, double z0, double h, double eps_r, double freq)
{
	double k, z2, z3, z2t, z3t, r;
	double w2, w3, w2t, w3t;
	double lambda2, lambda3, lambda2t, lambda3t;
	double tmp_eps_eff;

	printf(
		STR_PREFIX"Wilkinson:\n\n"
		"            +------+----[Z2t]-----o  \n"
		"           /       |                 \n"
		"         [Z2]      |           Port 2\n"
		"         /         ^                 \n"
		"    o---+          R                 \n"
		"Port 1   \\         v                 \n"
		"         [Z3]      |           Port 3\n"
		"           \\       |                 \n"
		"            +------+----[Z3t]-----o  \n"
		"\n\n"
			);

	k = sqrt(k2);

	z2 = z0 * sqrt(1/(k*k*k) + 1/k);
	z3 = k2 * z2;
	r = z0 * (k + 1/k);
	z2t = z0 / sqrt(k);
	z3t = z0 * sqrt(k);

	w2 = z2w(h, eps_r, z2);
	w3 = z2w(h, eps_r, z3);
	w2t = z2w(h, eps_r, z2t);
	w3t = z2w(h, eps_r, z3t);

	w2z(h, eps_r, w2, &tmp_eps_eff);
	lambda2 = calc_lambda(freq, tmp_eps_eff) / 1e6;
	w2z(h, eps_r, w3, &tmp_eps_eff);
	lambda3 = calc_lambda(freq, tmp_eps_eff) / 1e6;
	w2z(h, eps_r, w2t, &tmp_eps_eff);
	lambda2t = calc_lambda(freq, tmp_eps_eff) / 1e6;
	w2z(h, eps_r, w3t, &tmp_eps_eff);
	lambda3t = calc_lambda(freq, tmp_eps_eff) / 1e6;

	printf(
		STR_PREFIX"Z2 = %lf Ohm ,\tW = %lf mm ,\t Length = %lf mm\n"
		STR_PREFIX"Z3 = %lf Ohm ,\tW = %lf mm ,\t Length = %lf mm\n"
		STR_PREFIX"Lumped Resistor: R = %lf Ohm\n"
		STR_PREFIX"Z2t = %lf Ohm ,\tW = %lf mm ,\t Length = %lf mm\n"
		STR_PREFIX"Z3t = %lf Ohm ,\tW = %lf mm ,\t Length = %lf mm\n",
		z2, w2, lambda2/4,
		z3, w3, lambda3/4,
		r,
		z2t, w2t, (k == 1) ? 0 : lambda2t/4,
		z3t, w3t, (k == 1) ? 0 : lambda3t/4
	);
}

void novel_loop(double k2, double z0, double h, double eps_r, double freq)
{
	double z;
	double l_cos, w;
	double lambda;
	double tmp_eps_eff;

	printf(
		STR_PREFIX"A novel loop:\n\n"
		"           +-----------[Z,L2]-------------+         \n"
		"           |                              |   Port 2\n"
		"   o-------+                              +------o  \n"
		"Port 1     |                              |         \n"
		"           |                             <R>        \n"
		"           |                              |         \n"
		"           +----[Z,L3]----+----[Z,L3']----+         \n"
		"                          |                         \n"
		"                          |                         \n"
		"                          |                         \n"
		"                Port 3    o                         \n"
		"\n\n"
		);

	z = z0 * sqrt(2);

	w = z2w(h, eps_r, z);

	w2z(h, eps_r, w, &tmp_eps_eff);
	lambda = calc_lambda(freq, tmp_eps_eff) / 1e6;
	l_cos = 0.5*lambda * acos(1/sqrt(k2)) / M_PI;

	printf(
		STR_PREFIX"R = 100 Ohm ,\t Z = %lf Ohm ,\t W = %lf mm\n"
		STR_PREFIX"L2 = %lfmm ,\t L3 = %lf mm ,\t L3' = %lfmm\n"
		, z, w, lambda/4 + l_cos, lambda/4, l_cos);
}
/*****************************************************************************/

static double read_eps_r()
{
	double eps_r;

	printf("Relative epsilon = ");
	scanf("%lf", &eps_r);

	return eps_r;
}

static double read_height()
{
	double h;

	printf("Height of the substrat (in mm) = ");
	scanf("%lf", &h);

	return h;
}

static double read_k2()
{
	double k2;

	printf(
		"       P2      |S21|\n"
		"k^2 = ---- = (-------) ^ 2 = \n"
		"       P3      |S31|           "
		);
	scanf("%lf", &k2);

	return k2;
}

static double read_z0()
{
	double z0;

	printf("Z0 (in Ohm) = ");
	scanf("%lf", &z0);

	return z0;
}

static double read_freq()
{
	double freq;

	printf("Frequency (in GHz) = ");
	scanf("%lf", &freq);

	return freq;
}

/*****************************************************************************/

int main()
{
	double freq, h, eps_r, eps_r_eff, w0, z0, k2, lambda;

	freq = read_freq();
	eps_r = read_eps_r();
	h = read_height();
	z0 = read_z0();
	k2 = read_k2();

	BAR();
	w0 = z2w(h, eps_r, z0);

	printf (STR_PREFIX
		"Numerical tolerance of calculation for Z0 is %lf%%\n",
		100.0 * (w2z(h, eps_r, w0, &eps_r_eff) - z0) / z0);
	printf(STR_PREFIX"Width for %lf Ohm = %lf mm\n", z0, w0);

	lambda = calc_lambda(freq, eps_r_eff) / 1e6;

	printf(STR_PREFIX"Lambda for %lf Ohm at %lf GHz is %lf mm\n",
		       	z0, freq, lambda);
	printf(STR_PREFIX"Lambda / 2 is %lf mm\n", lambda/2);
	printf(STR_PREFIX"Lambda / 4 is %lf mm\n", lambda/4);
	BAR();

	branchline(k2, z0, h, eps_r, freq);
	BAR();

	wilkinson(k2, z0, h, eps_r, freq);
	BAR();

	novel_loop(k2, z0, h, eps_r, freq);
	BAR();

	return 0;
}
