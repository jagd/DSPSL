#include <stdio.h>
#include <stdlib.h>

#include "md.h"

#define error printf

#ifdef MD_ENABLE_DEBUG
void md_dump(struct MD *m)
{
	int row, col, i;

	i = 0;
	for (row = 0; row < m->rows; ++row) {
		for (col = 0; col < m->cols; ++col) {
			// printf("%12.3le", m->buf[i++]);
			printf("%30.20le", m->buf[i++]);
		}
		putchar('\n');
	}
}
#endif

struct MD* md_init(int rows, int cols)
{
	struct MD *m;

	if (rows <= 0 || cols <= 0) {
		error("md_init():"
			"number of rows and columns must be positive\n");
	}

	m = (struct MD*)malloc(sizeof(struct MD));
	m->rows = rows;
	m->cols = cols;
	m->buf = (double*)malloc(sizeof(double)*rows*cols);

	if ((m == 0) || (m->buf == 0)) {
		error("md_init(): not enough memory\n");
	}

	return m;
}


void md_fill(struct MD *m, double val)
{
	int i, t;

	t = m->cols * m->rows;

	for (i = 0; i < t; ++i) {
		m->buf[i] = val;
	}
}


void md_free(struct MD *m)
{
	if (m) {
		if (m->buf) {
			free(m->buf);
		}
		free(m);
	}
}


void md_set(struct MD *m, int row, int col, double val)
{
#ifdef MD_ENABLE_RANGE_CHECKING
	if ((row < 0 || row >= m->rows)
		|| (col < 0 || col >= m->cols)) {
		error("md_set(): index out of range\n");
		return;
	}
#endif
	m->buf[row*m->cols + col] = val;
}


double md_get(struct MD *m, int row, int col)
{
#ifdef MD_ENABLE_RANGE_CHECKING
	if ((row < 0 || row >= m->rows)
		|| (col < 0 || col >= m->cols)) {
		error("md_set(): index out of range\n");
		return 0;
	}
#endif
	return m->buf[row*m->cols + col];
}


void md_inverse_direct(struct MD *m)
{
	int *row_exchange;
	int row, col;

#ifdef MD_ENABLE_RANGE_CHECKING
	if (m->rows != m->cols) {
		error("md_inverse_direct(): Matrix must be square");
		return;
	}
#endif

	row_exchange = (int*)malloc(m->rows * sizeof(int));

	for (row = 0; row < m->rows; ++row) {
		/************* 1 *************/
		/* relative column pivoting */
		/*
		   not from the first column,
		   just as the elimination,
		   treating the first columns were 0
		*/

		double factor;

		double pivot_val = 0; /* the pivot value */

		int row2;
		double col_max = 0;

		row_exchange[row] = row; /* init */

		for (row2 = row; row2 < m->rows; ++row2) {
			int sum;
			double curr;
			double unscale_val;

			unscale_val = m->buf[row2*m->cols + row];
			if (unscale_val == 0) {
				continue;
			}

			sum = 0;
			for (col = row; col < m->cols; ++col) {
				int val = m->buf[row2*m->cols + col];
				sum += val >= 0 ? val : -val; /* inline abs */
			}

			curr = unscale_val / sum;
			curr = curr >= 0 ? curr : -curr; /* inline abs*/

			if (curr > col_max) {
				pivot_val = unscale_val;
				col_max = curr;
				row_exchange[row] = row2;
			}
		}

		/* do rows exchange, move the pivoting row to the next */

		if (row != row_exchange[row]) {
			/* for 2 row operations, to save the multiplication */
			int base1, base2;

			base1 = row*m->cols;
			base2 = row_exchange[row]*m->cols;
			for (col = 0; col < m->cols; ++col) {
				double tmp;

				tmp = m->buf[base1 + col];
				m->buf[base1 + col] = m->buf[base2 + col];
				m->buf[base2 + col] = tmp;
			}
		}

#ifdef MD_ENABLE_DEBUG
		if (pivot_val == 0) {
			error("md_inverse_direct(): Matrix is singular\n");
			return;
		}
#endif
		/************* 2 *************/
		/* the elimination */

		/* first, calculate the exchanged equation */
		factor = 1/pivot_val;
		m->buf[row*m->cols + row] = factor;
		factor = -factor;
		for (col = 0; col < row; ++col) {
			m->buf[row*m->cols + col] *= factor;
		}
		for (col = row + 1; col < m->cols; ++col) {
			m->buf[row*m->cols + col] *= factor;
		}

		/* second, the other equation */
		/* except col == row */
		for (col = 0; col < m->cols; ++col) {
			if (col == row) {
				continue;
			}


			for (row2 = 0; row2 < row; ++row2) {
				factor = m->buf[row2*m->cols + row]
						* m->buf[row*m->cols + col];
				m->buf[row2*m->cols + col] += factor;
			}
			for (row2 = row + 1; row2 < m->rows; ++row2) {
				factor = m->buf[row2*m->cols + row]
						* m->buf[row*m->cols + col];
				m->buf[row2*m->cols + col] += factor;
			}
		}
		/* for col == row */
		factor = 1/pivot_val;
		for (row2 = 0; row2 < row; ++row2) {
			m->buf[row2*m->cols + row] *= factor;
		}
		for (row2 = row + 1; row2 < m->rows; ++row2) {
			m->buf[row2*m->cols + row] *= factor;
		}

	}

	/* column exchanging (not row !) */
	for (col = m->rows - 2; col >= 0; --col) {
		int col2 = row_exchange[col];
		if (col == col2) {
			continue;
		}

		for (row = 0; row < m->rows; ++row) {
			double tmp;
			tmp = m->buf[row*m->cols + col];
			m->buf[row*m->cols + col] =
				m->buf[row*m->cols + col2];
			m->buf[row*m->cols + col2] = tmp;
		}
	}

	free(row_exchange);
}

struct MD* md_mul(struct MD *a, struct MD *b)
{
	struct MD *c;

	int row, col;

#ifdef MD_ENABLE_RANGE_CHECKING
	if (a->cols != b->rows) {
		error("md_mul(): Uncorrect matrix dimension\n");
	}
#endif

	c = md_init(a->rows, b->cols);

	for (row = 0; row < c->rows; ++row) {
		for (col = 0; col < c->cols; ++col) {
			int i;
			c->buf[row*c->cols + col] = 0;
			for (i = 0; i < a->cols; ++i) {
				c->buf[row*c->cols + col] +=
					a->buf[row*a->cols + i]
					* b->buf[i*b->cols + col];
			}
		}
	}

	return c;
}
