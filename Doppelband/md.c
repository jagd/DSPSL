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

		double factor;

		double pivot_val; /* the pivot value */

		int row2;
		double col_max = 0;

		row_exchange[row] = row; /* init */

		for (row2 = row; row2 < m->rows; ++row2) {
			int index; /* to save some multiplications */
			int sum;
			double curr;
			double unscale_val;

			/* row2 <+> the outer column index */
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
				pivot_val = unscale_val;
				col_max = curr;
				row_exchange[row] = row2;
			}
		}

		/* do rows exchange, move the pivoting row to the next */

		if (row != row_exchange[row]) {
			/* for 2 row operations, to save the multiplication */
			int index1, index2;

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
		}

printf("after exchange for row %d:\n", row);
md_dump(m);

		/************* 2 *************/
		/* the elimination */

printf("pivot_val = %lf:\n", pivot_val);

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

printf("between elimination for row %d:\n", row);
md_dump(m);
		/* second, the other equation */
		/* except col == row */
		for (col = 0; col < m->cols; ++col) {
			if (col == row) {
				continue;
			}


			for (row2 = 0; row2 < row; ++row2) {
				factor = m->buf[row2*m->cols + row] * m->buf[row*m->cols + col];
				m->buf[row2*m->cols + col] += factor;
			}
			for (row2 = row + 1; row2 < m->rows; ++row2) {
				factor = m->buf[row2*m->cols + row] * m->buf[row*m->cols + col];
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
printf("after elimination for row %d:\n", row);
md_dump(m);

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

int main()
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

	return 0;
}
