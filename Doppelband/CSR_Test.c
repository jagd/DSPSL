#include <stdio.h>
#include <stdlib.h>

#include "CSR.h"

/******************************************/
/*              Unit Test                 */
/******************************************/
static void dump_readable(struct CSR* m)
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
			printf("%10.2f ", CSR_get(m, i, j));
		}
		puts("");
	}

	puts("\n");
}

static void dump_raw(struct CSR* m)
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

int main()
{
	struct CSR* m;
	struct CSR* a;

	// dump_raw(NULL);

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

	printf("\n\nAddition:\n");
	CSR_add_online(m, m); /* overwrite */
	dump_readable(m);

	printf("\n\nAddition with:\n");
	a = CSR_init(3, 4, 100);
	CSR_set(a, 0, 0, 1);
	CSR_set(a, 1, 1, 2);
	CSR_set(a, 2, 2, 3);
	CSR_set(a, 2, 3, 4);
	CSR_set(a, 0, 3, 5);
	dump_readable(a);
	printf("\nResult :\n");
	CSR_add_online(a, m); /* overwrite */
	dump_readable(a);

	CSR_free(a);
	CSR_free(m);
	return 0;
}
