.PHONY: all clean

CFLAGS:=-O3 -Wall $(CFLAGS)

all: md_test mom

md_test: md_test.c md.c
	$(CC) -D MD_ENABLE_DEBUG  $(CFLAGS) $^ -o $@

mom: mom.c mom_mesh.c md.c
	$(CC) -D MOM_MESH_ENABLE_DEBUG -lm $(CFLAGS) $^ -o $@

clean:
	-rm *.o
	-rm mom md_test
