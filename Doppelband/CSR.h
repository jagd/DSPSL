#ifndef __CSR_H__
#define __CSR_H__

#define CSR_CHECK_DIMENSION 1

#define INCREASE_SIZE 1000

/* Should not be size_t , because a signed type is necessary */
/* In some case, say row or column > 32767 , 32 bit is not enough*/
typedef int INDEX_TYPE;


struct CSR {
	float *buf;
	INDEX_TYPE *row_index; /* size == row numbers as index */
	INDEX_TYPE *col_index; /* as many as buf[] */
	INDEX_TYPE allocated, size; /* size counter for buf[] */
	INDEX_TYPE rows, cols;
};

struct CSR* CSR_init(INDEX_TYPE rows, INDEX_TYPE cols, INDEX_TYPE pre_allocate);
void CSR_free(struct CSR* p);
void CSR_extend(struct CSR* m, INDEX_TYPE increase);
float CSR_get(struct CSR* m, INDEX_TYPE row, INDEX_TYPE col);
void CSR_set(struct CSR* m, INDEX_TYPE row, INDEX_TYPE col, float val);
void CSR_set_add(struct CSR* m, INDEX_TYPE row, INDEX_TYPE col, float toadd);
void CSR_add_online(struct CSR *a, struct CSR *b);

#endif
