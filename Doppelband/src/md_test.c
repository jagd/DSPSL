#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "md.h"

void inv_big_random()
{
	int size = 500;
	struct MD *a;
	int i, j;

	a = md_init(size, size);

	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			md_set(a, i, j, random()%2000 - 1000);
		}
	}

	md_inverse_direct(a);

	md_free(a);
}

void inv_small_random()
{
	int size = 50;
	struct MD *a;
	int i, j;

	a = md_init(size, size);

	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			md_set(a, i, j, random()%2000 - 1000);
		}
	}

	md_inverse_direct(a);

	md_free(a);
}

void inv_tiny()
{
	struct MD *a;

	a = md_init(3, 3);
	md_fill(a, 0);
	md_set(a, 0, 0, 2);
	md_set(a, 0, 1, 7);
	md_set(a, 1, 0, 5);
	md_set(a, 1, 1, 3);

	md_set(a, 0, 2, -1);
	md_set(a, 1, 2, -5);
	md_set(a, 2, 0, 11);
	md_set(a, 2, 1, 12);
	md_set(a, 2, 2, 22);

	md_dump(a);

	md_inverse_direct(a);
	printf("\ninversed:\n");
	md_dump(a);

	md_free(a);
}

void mul_tiny()
{
	struct MD *a, *b, *c;

	a = md_init(2, 2);
	md_fill(a, 0);
	md_set(a, 0, 0, 1);
	md_set(a, 0, 1, 2);
	md_set(a, 1, 0, 3);
	md_set(a, 1, 1, 4);

	b = md_init(2, 1);
	md_fill(b, 0);
	md_set(b, 0, 0, 2);
	md_set(b, 1, 0, 3);


	printf("\nA =\n");
	md_dump(a);

	printf("\nB =\n");
	md_dump(b);

	printf("\nMultiplication result:\n");
	c = md_mul(a, b);
	md_dump(c);

	md_free(a);
	md_free(b);
	md_free(c);
}

void eye(int n)
{
	struct MD *a;

	a = md_eye(n);

	printf("\neye :\n");
	md_dump(a);

	md_free(a);
}

int main()
{
	eye(5);

	printf("=== inversion tiny ===\n");
	inv_tiny();
	printf("=== inversion small ===\n");
	inv_small_random();
	printf("=== inversion big ===\n");
	inv_big_random();

	printf("=== mul tiny ===\n");
	mul_tiny();

	return 0;
}
