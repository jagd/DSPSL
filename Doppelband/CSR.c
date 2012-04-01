#include <stdio.h>
#include <stdlib.h>

#include "CSR.h"

#define error printf

struct CSR* CSR_init(INDEX_TYPE rows, INDEX_TYPE cols, INDEX_TYPE pre_allocate)
{
	struct CSR* m;

	if (pre_allocate == 0 || rows == 0) {
		error("pre_allocate and row_index must bigger than 0\n");
		return NULL;
	}

	m = (struct CSR*) malloc(sizeof(struct CSR));

	m->rows = rows;
	m->cols = cols;

	m->size = 0;
	m->row_index = (INDEX_TYPE*) calloc(rows + 1, sizeof(INDEX_TYPE));

	m->allocated = pre_allocate;
	m->buf = (float*) malloc(pre_allocate * sizeof(float));
	m->col_index = (INDEX_TYPE*) malloc(pre_allocate * sizeof(INDEX_TYPE));

	return m;
}


void CSR_free(struct CSR* p)
{
	free(p->row_index);
	free(p->buf);
	free(p->col_index);
	free(p);
}


void CSR_extend(struct CSR* m, INDEX_TYPE increase)
{
	m->allocated += increase;

	m->buf = realloc(m->buf, (m->allocated) * sizeof(float));

	m->col_index = realloc(m->col_index,
			(m->allocated + 1) * sizeof(INDEX_TYPE));
}

float CSR_get(struct CSR* m, INDEX_TYPE row, INDEX_TYPE col)
{
	INDEX_TYPE left, right;
	/* index in the buffer, left <= the target < right */
	INDEX_TYPE mid = 0;

	left = m->row_index[row];
	right = m->row_index[row + 1];

	if (left == right) {
		return 0;
	}

	while (left < right) {
		/* left <= target < right */
		mid = (left + right) >> 1;
		if (m->col_index[mid] == col) {
			break;
		}
		if (m->col_index[mid] < col) {
			left = mid + 1;
		} else {
			right = mid;
		}
	}

	if (left == right) {
		return 0;
	} else {
		return m->buf[mid];
	}

}

void CSR_set(struct CSR* m, INDEX_TYPE row, INDEX_TYPE col, float val)
{
	INDEX_TYPE left, right;
	/* index in the buffer, left <= the target < right */
	INDEX_TYPE mid = 0;

	left = m->row_index[row];
	right = m->row_index[row + 1];

	if (left == right) {
		/* row not found */
		/* memmove() from string.h has overhead for a temporary area */
		INDEX_TYPE i;

		if (m->size == m->allocated) {
			CSR_extend(m, INCREASE_SIZE);
		}

		for (i = m->size; i >= left; --i) {
			m->buf[i] = m->buf[i-1];
			m->col_index[i] = m->col_index[i-1];
		}
		m->buf[left] = val;
		m->col_index[left] = col;

		/* this step exist only by row inserting */
		m->row_index[row] = left;

		for (i = row + 1; i <= m->rows; ++i) {
			m->row_index[i]++;
		}
		m->size++;
		return;
	}

	while (left < right) {
		/* left <= target < right */
		mid = (left + right) >> 1;
		if (m->col_index[mid] == col) {
			break;
		}
		if (m->col_index[mid] < col) {
			left = mid + 1;
		} else {
			right = mid;
		}
	}

	if (left == right) {
		/* col not found */
		INDEX_TYPE i;

		if (m->size == m->allocated) {
			CSR_extend(m, INCREASE_SIZE);
		}

		for (i = m->size; i >= left; --i) {
			m->buf[i] = m->buf[i-1];
			m->col_index[i] = m->col_index[i-1];
		}
		m->buf[left] = val;
		m->col_index[left] = col;

		for (i = row + 1; i <= m->rows; ++i) {
			m->row_index[i]++;
		}
		m->size++;
	} else {
		/* there was a element at the same place, just replace it */
		m->buf[mid] = val;
	}
}

void CSR_set_add(struct CSR* m, INDEX_TYPE row, INDEX_TYPE col, float toadd)
{
	INDEX_TYPE left, right;
	/* index in the buffer, left <= the target < right */
	INDEX_TYPE mid = 0;

	left = m->row_index[row];
	right = m->row_index[row + 1];

	if (left == right) {
		/* row not found */
		/* memmove() from string.h has overhead for a temporary area */
		INDEX_TYPE i;

		if (m->size == m->allocated) {
			CSR_extend(m, INCREASE_SIZE);
		}

		for (i = m->size; i >= left; --i) {
			m->buf[i] = m->buf[i-1];
			m->col_index[i] = m->col_index[i-1];
		}
		m->buf[left] = toadd;
		m->col_index[left] = col;

		/* this step exist only by row inserting */
		m->row_index[row] = left;

		for (i = row + 1; i <= m->rows; ++i) {
			m->row_index[i]++;
		}
		m->size++;
		return;
	}

	while (left < right) {
		/* left <= target < right */
		mid = (left + right) >> 1;
		if (m->col_index[mid] == col) {
			break;
		}
		if (m->col_index[mid] < col) {
			left = mid + 1;
		} else {
			right = mid;
		}
	}

	if (left == right) {
		/* col not found */
		INDEX_TYPE i;

		if (m->size == m->allocated) {
			CSR_extend(m, INCREASE_SIZE);
		}

		for (i = m->size; i >= left; --i) {
			m->buf[i] = m->buf[i-1];
			m->col_index[i] = m->col_index[i-1];
		}
		m->buf[left] = toadd;
		m->col_index[left] = col;

		for (i = row + 1; i <= m->rows; ++i) {
			m->row_index[i]++;
		}
		m->size++;
	} else {
		/* there was a element at the same place, just replace it */
		m->buf[mid] += toadd;
	}
}

/* A := A + B */
void CSR_add_online(struct CSR *a, struct CSR *b)
{
	INDEX_TYPE row;

#ifdef CSR_CHECK_DIMENSION
	if ((a->rows != b->rows)
		|| (a->cols != b->cols))
	{
		error("CSR_add_online(): Matrixs must have the same dimension");
	}
#endif
	/* for all in B, add it to A */

	for (row = 0; row < b->rows; ++row) {
		INDEX_TYPE i;
		INDEX_TYPE right = b->row_index[row + 1];
		/* row exist */
		for (i = b->row_index[row]; i < right; ++i) {
			CSR_set_add(a, row, b->col_index[i], b->buf[i]);
		}
	}
}
