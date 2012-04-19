.PHONY: all clean

CFLAGS:=-O3 -Wall $(CFLAGS)

all: md_test mom_test

md_test: md_test.c md.c
	$(CC) -D MD_ENABLE_DEBUG  $(CFLAGS) $^ -o $@

mom_test: mom_test.c mom.c mom_mesh.c md.c
	$(CC) -D MD_ENABLE_DEBUG -D MOM_MESH_ENABLE_DEBUG -lm $(CFLAGS) $^ -o $@

clean:
	-rm *.o
	-rm mom_test md_test
