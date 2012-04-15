#ifndef __MD_H__
#define __MD_H__

/* #define MD_ENABLE_DEBUG 1 */

/* Matrix of Double */
struct MD {
	int rows, cols;
	double *buf;
};

struct MD* md_init(int rows, int cols);
void md_fill(struct MD *m, double val);
void md_free(struct MD *m);

/* not inline: */
void md_set(struct MD *m, int row, int col, double val);
double md_get(struct MD *m, int row, int col);

void md_inverse_direct(struct MD *m);

#ifdef MD_ENABLE_DEBUG
#define MD_ENABLE_RANGE_CHECKING 1
void md_dump(struct MD *m);
#endif


#endif
