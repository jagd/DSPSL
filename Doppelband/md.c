#include <stdio.h>
#include <stdlib.h>

#define error printf
#define INLINE inline

#define ENABLE_RANGE_CHECKING 1

/* Matrix of Double */
struct MD {
	int rows, cols;
	double *buf;
};

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

void INLINE md_set(struct MD *m, int row, int col, double val)
{
#ifdef ENABLE_RANGE_CHECKING
	if ((row < 0 || row >= m->rows)
		|| (col < 0 || col >= m->cols)) {
		error("md_set(): index out of range\n");
	}
#endif
	m->buf[row*m->cols + col] = val;
}

double INLINE md_get(struct MD *m, int row, int col)
{
#ifdef ENABLE_RANGE_CHECKING
	if ((row < 0 || row >= m->rows)
		|| (col < 0 || col >= m->cols)) {
		error("md_set(): index out of range\n");
	}
#endif
	return m->buf[row*m->cols + col];
}

void md_dump(struct MD *m)
{
	int row, col, i;

	i = 0;
	for (row = 0; row < m->rows; ++row) {
		for (col = 0; col < m->cols; ++col) {
			printf("%12.3lf", m->buf[i++]);
		}
		putchar('\n');
	}
}

void md_inverse_direct(struct MD *m)
{
	int *row_exchange;
	int row, col;

	row_exchange = (int*)malloc(m->rows * sizeof(int));

	for (row = 0; row < m->rows; ++row) {
		/************* 1 *************/
		/* relative column pivoting */
		/*
		   not from the first column,
		   just as the elimination,
		   treating the first columns were 0
		*/

		int norm; /* the pivot value */

		/* for 2 row operations, to save the multiplication */
		int index1, index2;

		int row2;
		double col_max = 0;

		row_exchange[row] = row; /* init */

		for (row2 = row; row2 < m->rows; ++row2) {
			int index; /* to save some multiplications */
			int sum;
			double curr;
			double unscale_val;

			/* row2 <+> the outter column index */
			index = row2*m->cols + row;

			/* to save the `sum` calculation */
			unscale_val = m->buf[index];
			if (unscale_val == 0) {
				continue;
			}

			sum = 0;
			for (col = row; col < m->cols; ++col) {
				int val = m->buf[index];
				sum += val >= 0 ? val : -val; /* inline abs */
				++index;
			}

			curr = unscale_val / sum;
			curr = curr >= 0 ? curr : -curr; /* inline abs*/

			if (curr > col_max) {
				norm = unscale_val;
				col_max = curr;
				row_exchange[row] = row2;
			}
		}

		/* do rows exchange, move the pivoting row to the next */

		/* to save multiplications */
		index1 = row*m->cols;
		index2 = row_exchange[row]*m->cols;
		for (col = 0; col < m->cols; ++col) {
			double tmp;

			tmp = m->buf[index1];
			m->buf[index1] = m->buf[index2];
			m->buf[index2] = tmp;

			++index1;
			++index2;
		}

		/************* 2 *************/
		/* the elimination */
		// TODO: use norm
		for (j = 0; j < m->cols; ++j) {
		}

	}

	/* column exchanging (not row !) */
	// TODO: 

	free(row_exchange);
}

int main()
{

	struct MD *a;

	a = md_init(3, 3);
	md_fill(a, 0);
	md_set(a, 0, 0, 1);
	md_set(a, 1, 1, 2);
	md_set(a, 2, 2, 1);

	md_dump(a);

	md_free(a);

	return 0;
}
