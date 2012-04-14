#include <stdio.h>
#include <stdlib.h>

#include "CSR.h"

#define error printf
#define EPSILON 1e-5

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
	if (p == NULL) {
		return;
	}

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

#ifdef CSR_CHECK_DIMENSION
	if (row < 0 || row >= m->rows || col < 0 || col >= m->cols) {
		error("CSR_get(): invalid row or col");
	}
#endif

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

#ifdef CSR_CHECK_DIMENSION
	if (row < 0 || row >= m->rows || col < 0 || col >= m->cols) {
		error("CSR_set(): invalid row or col");
	}
#endif

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

#ifdef CSR_CHECK_DIMENSION
	if (row < 0 || row >= m->rows || col < 0 || col >= m->cols) {
		error("CSR_set_add(): invalid row or col");
	}
#endif

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


/*
	C = AB
	Sparse A can get the best performance.
	WARNING: the insert of C can be very inefficient!
*/
struct CSR* CSR_mul(struct CSR *a, struct CSR* b)
{
	struct CSR* c;
	INDEX_TYPE row;

#ifdef CSR_CHECK_DIMENSION
	if (a->cols != b->rows) {
		error("CSR_mul(A,B): dimension of A and B not correct.\n");
		return NULL;
	}
#endif

	c = CSR_init( (a->rows < b->rows) ? a->rows : b->rows,
			(a->cols < b->cols) ? a->cols : b->cols,
			(a->size > b->size)?a->size:b->size );

	/* forall rows in a --> rows in c with multiplication col in b */
	for (row = 0; row < a->rows; ++row) {
		INDEX_TYPE col_c;
		INDEX_TYPE left, right;
		left = a->row_index[row];
		right = a->row_index[row + 1];
		if (left >=  right) {
			/* row does not exist => c[ _ , row ] = 0 */
			continue;
		}
		/* row exist */
		for (col_c = 0; col_c < b->cols; ++col_c) {
			INDEX_TYPE i_a;
			float v = 0;
			for (i_a= left; i_a < right; ++i_a) {
				v += a->buf[i_a]
					* CSR_get(b, a->col_index[i_a], col_c);
			}
			if (v != 0) {
				CSR_set(c, row, col_c, v);
			}
		}
	}

	return c;
}


/*
	calculate `x` from

		Ax = b
	where `b` is a (vertical) vector

	A = L + D + U
==>
	x = D^(-1) * (b - (L + U) * x)
*/
/* !!!FIXME!!! */
struct CSR* CSR_Jacobi(struct CSR *a, struct CSR *b)
{
	struct CSR *x, *x_old;
	INDEX_TYPE row;
	float eps = 1e37;

#ifdef CSR_CHECK_DIMENSION
	if (b->cols != 1) {
		error("Jacobi(A, B): B should have only one column.\n");
	}
	if ((a->cols != b->rows) || (a->rows != a->cols)) {
		error("Jacobi(A, B): dimension of A and B not correct.\n");
	}
#endif

	x = CSR_init(b->rows, 1, b->rows);
	x_old = CSR_init(b->rows, 1, b->rows);
	/* dirty filling: init all elements to 0 */
	for (row = 0; row < x->rows; ++row) {
		x->buf[row] = 0; /* start value */
		x->row_index[row] = row;
		x->col_index[row] = 0;
		x_old->buf[row] = 0; /* start value */
		x_old->row_index[row] = row;
		x_old->col_index[row] = 0;
	}
	x->size = b->rows;
	x->row_index[row] = row; /* last fake row_index */
	x_old->size = b->rows;
	x_old->row_index[row] = row; /* last fake row_index */

	while (eps > EPSILON) {
		for (row = 0; row < x->rows; ++row) {
			x_old->buf[row] = x->buf[row];
		}

		/* forall rows in A */

		for (row = 0; row < x->rows; ++row) {
			INDEX_TYPE left, right, i;
			float val;

			left = a->row_index[row];
			right = a->row_index[row + 1];
			if (left >= right) {
				/* row does not exist */
				continue;
			}
			val = 0;
			for (i = left; i < right; ++i) {
				if (row == a->col_index[i]) {
					continue;
				}
				val -= a->buf[i] * x->buf[a->col_index[i]];
			}
			x->buf[row] = val;
		}

		CSR_add_online(x, b);

		for (row = 0; row < x->rows; ++row) {
			x->buf[row] /= CSR_get(a, row, row);
		}

		eps = 0;
		for (row = 0; row < x->rows; ++row) {
			float delta;
			delta = x->buf[row] - x_old->buf[row];
			eps += (delta > 0) ? delta : -delta;
		}
	}

	CSR_free(x_old);
	return x;
}

/*
	calculate `x` from

		Ax = b
	where `b` is a (vertical) vector

	A = L + D + U
==>
	x = D^(-1) * (b - (L + U) * x)
*/
struct CSR* CSR_SOR(struct CSR *a, struct CSR *b)
{
	struct CSR *x, *x_old;
	INDEX_TYPE row;
	float eps = 1e37;

#ifdef CSR_CHECK_DIMENSION
	if (b->cols != 1) {
		error("Jacobi(A, B): B should have only one column.\n");
	}
	if ((a->cols != b->rows) || (a->rows != a->cols)) {
		error("Jacobi(A, B): dimension of A and B not correct.\n");
	}
#endif

	x = CSR_init(b->rows, 1, b->rows);
	x_old = CSR_init(b->rows, 1, b->rows);
	/* dirty filling: init all elements to 0 */
	for (row = 0; row < x->rows; ++row) {
		x->buf[row] = 0; /* start value */
		x->row_index[row] = row;
		x->col_index[row] = 0;
		x_old->buf[row] = 0; /* start value */
		x_old->row_index[row] = row;
		x_old->col_index[row] = 0;
	}
	x->size = b->rows;
	x->row_index[row] = row; /* last fake row_index */
	x_old->size = b->rows;
	x_old->row_index[row] = row; /* last fake row_index */

	while (eps > EPSILON) {
		for (row = 0; row < x->rows; ++row) {
			x_old->buf[row] = x->buf[row];
		}

		/* forall rows in A */

		for (row = 0; row < x->rows; ++row) {
			INDEX_TYPE left, right, i;
			float val;

			left = a->row_index[row];
			right = a->row_index[row + 1];
			if (left >= right) {
				/* row does not exist */
				continue;
			}
			val = 0;
			for (i = left; i < right; ++i) {
				if (row == a->col_index[i]) {
					continue;
				}
				if (row < a->col_index[i]) {
					val -= a->buf[i] * x->buf[a->col_index[i]];
				} else {
					val -= a->buf[i] * x_old->buf[a->col_index[i]];
				}
			}
			x->buf[row] = val;
		}

		CSR_add_online(x, b);

		for (row = 0; row < x->rows; ++row) {
			x->buf[row] *= CSR_SOR_OMEGA;
			x->buf[row] /= CSR_get(a, row, row);
			x->buf[row] += (1-CSR_SOR_OMEGA) * x_old->buf[row];
		}

		eps = 0;
		for (row = 0; row < x->rows; ++row) {
			float delta;
			delta = x->buf[row] - x_old->buf[row];
			eps += (delta > 0) ? delta : -delta;
		}
	}

	CSR_free(x_old);
	return x;
}
