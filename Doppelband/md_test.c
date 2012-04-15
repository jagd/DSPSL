#include <stdio.h>
#include <stdlib.h>

#include "md.h"

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
