#include <stdio.h>
#include <stdlib.h>

#include "CSR.h"

/******************************************/
/*              Unit Test                 */
/******************************************/
void dump_readable(struct CSR* m)
{
	INDEX_TYPE i, j;

	if (m == NULL) {
		printf("dump_readable(): NULL point for the array.\n");
		return;
	}

	printf("\n========== Readable Dump ==========\n");
	printf("Elements = %d\n", m->size);
	printf("Dimension = %d x %d\n", m->rows, m->cols);

	printf("\nValues:\n");

	for (i = 0; i < m->rows; ++i) {
		for (j = 0; j < m->cols; ++j) {
			printf("%10.4f ", CSR_get(m, i, j));
		}
		puts("");
	}

	puts("\n");
}

void dump_raw(struct CSR* m)
{
	INDEX_TYPE i;

	if (m == NULL) {
		printf("dump_raw(): NULL point for the array.\n");
		return;
	}

	printf("\n========== Raw Dump ==========\n");
	printf("Elements = %d\n", m->size);
	printf("Dimension = %d x %d\n", m->rows, m->cols);

	printf("\nValues:\n");
	for (i = 0; i < m->size; ++i) {
		printf("%f\t", m->buf[i]);
	}

	printf("\nRow Index Vector:\n");
	for (i = 0; i < m->rows; ++i) {
		printf("%d\t", m->row_index[i]);
	}
	printf("|| %d", m->row_index[i]);  /* the last flag element */

	printf("\nColumn Vector:\n");
	for (i = 0; i < m->size; ++i) {
		printf("%d\t", m->col_index[i]);
	}

	putchar('\n');
	putchar('\n');
}

void test_create()
{
	struct CSR* m;

	m = CSR_init(3, 4, 1);

	CSR_set(m, 0, 0, 100);
	CSR_set(m, 0, 1, 2);
	CSR_set(m, 1, 1, 3);
	CSR_set(m, 1, 2, 9);
	CSR_set(m, 2, 1, 1);
	CSR_set(m, 2, 2, 4);

	CSR_set(m, 0, 0, 1); /* overwrite */

	dump_raw(m);
	dump_readable(m);

	CSR_free(m);
}

void test_add()
{
	struct CSR *a, *b;

	a = CSR_init(3, 4, 100);
	CSR_set(a, 0, 0, 1);
	CSR_set(a, 1, 1, 2);
	CSR_set(a, 2, 2, 3);
	CSR_set(a, 2, 3, 4);
	CSR_set(a, 0, 3, 5);
	dump_readable(a);

	b = CSR_init(3, 4, 1);
	CSR_set(b, 0, 0, 1);
	CSR_set(b, 0, 1, 2);
	CSR_set(b, 1, 1, 3);
	CSR_set(b, 1, 2, 9);
	CSR_set(b, 2, 1, 1);
	CSR_set(b, 2, 2, 4);
	dump_readable(a);

	CSR_add_online(a, b); /* overwrite */
	dump_readable(a);

	CSR_free(a);
	CSR_free(b);
}

void test_mul()
{
	struct CSR* a;
	struct CSR* b;
	struct CSR* c;

	a = CSR_init(3, 4, 100);
	CSR_set(a, 0, 0, 1);
	CSR_set(a, 1, 1, 2);
	CSR_set(a, 2, 2, 3);
	CSR_set(a, 2, 3, 4);
	CSR_set(a, 0, 3, 5);
	dump_readable(a);

	b = CSR_init(4, 3, 100);
	CSR_set(b, 3, 2, 10);
	CSR_set(b, 0, 0, 3);
	CSR_set(b, 1, 2, 4);
	CSR_set(b, 1, 1, 2);
	CSR_set(b, 3, 0, 5);
	CSR_set(b, 2, 1, 4);
	dump_readable(b);

	c = CSR_mul(a, b);
	dump_readable(c);

	CSR_free(a);
	CSR_free(b);
	CSR_free(c);
}

void test_mul_2()
{
	struct CSR* a;
	struct CSR* b;
	struct CSR* c;

	a = CSR_init(3, 4, 100);
	CSR_set(a, 0, 0, 1);
	CSR_set(a, 1, 1, 2);
	CSR_set(a, 2, 2, 3);
	CSR_set(a, 2, 3, 4);
	CSR_set(a, 0, 3, 5);
	dump_readable(a);

	b = CSR_init(4, 1, 100);
	CSR_set(b, 3, 0, 10);
	CSR_set(b, 2, 0, 4);
	CSR_set(b, 0, 0, 3);
	CSR_set(b, 1, 0, 4);
	dump_readable(b);

	c = CSR_mul(a, b);
	dump_readable(c);

	CSR_free(a);
	CSR_free(b);
	CSR_free(c);
}


void test_sor()
{
	struct CSR* a;
	struct CSR* b;
	struct CSR* x;

	a = CSR_init(4, 4, 10);
	CSR_set(a, 0, 0, 100);
	CSR_set(a, 0, 1, 2);
	CSR_set(a, 0, 2, 5);
	CSR_set(a, 0, 3, 4);
	CSR_set(a, 1, 0, 4);
	CSR_set(a, 1, 1, 470);
	CSR_set(a, 1, 2, 2);
	CSR_set(a, 1, 3, 9);
	CSR_set(a, 2, 0, 1);
	CSR_set(a, 2, 1, 8);
	CSR_set(a, 2, 2, 300);
	CSR_set(a, 2, 3, 1);
	CSR_set(a, 3, 0, 6);
	CSR_set(a, 3, 1, 3);
	CSR_set(a, 3, 2, 2);
	CSR_set(a, 3, 3, 500);
	dump_readable(a);

	b = CSR_init(4, 1, 100);
	CSR_set(b, 3, 0, 10);
	CSR_set(b, 2, 0, 7);
	CSR_set(b, 0, 0, 9);
	CSR_set(b, 1, 0, 4);
	dump_readable(b);

	x = CSR_SOR(a, b);
	dump_readable(x);

	CSR_free(a);
	CSR_free(b);
	CSR_free(x);
}

void test_sor_big(int n)
{
	struct CSR* a;
	struct CSR* b;
	struct CSR* x;
	int i,j;

	a = CSR_init(n, n, n*n);

	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			CSR_set(a, i, j, (rand() % 100000) / 10.0);
		}
		CSR_set(a, i, i, 100.0*CSR_get(a, i, i));
	}
	// dump_readable(a);

	b = CSR_init(n, 1, n);
	for (i = 0; i < n; ++i) {
		CSR_set(b, i, 0, 1000.0*(rand() % 100));
	}
	// dump_readable(b);

	x = CSR_SOR(a, b);
	dump_readable(x);

	CSR_free(a);
	CSR_free(b);
	CSR_free(x);
}

int main()
{
	srand(1000);
	test_sor_big(50);
	return 0;
}
