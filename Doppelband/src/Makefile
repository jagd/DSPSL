.PHONY: tags all ui clean test

CFLAGS:=-O3 -Wall $(CFLAGS)

all: mom.o mom_mesh.o md.o test ui

mom.o: mom.c
	$(CC) -c $(CFLAGS) $^ -o $@

mom_mesh.o: mom_mesh.c
	$(CC) -c $(CFLAGS) $^ -o $@

md.o: md.c
	$(CC) -c $(CFLAGS) $^ -o $@

test: md_test mom_test

md_test: md_test.c md.c
	$(CC) -D MD_ENABLE_DEBUG  $(CFLAGS) $^ -o $@

mom_test: mom_test.c mom.c mom_mesh.c md.c
	$(CC) -D MD_ENABLE_DEBUG -D MOM_MESH_ENABLE_DEBUG -lm $(CFLAGS) $^ -o $@

clean:
	-rm tags
	-rm *.o
	-rm mom_test md_test
	-rm *.plot
	-$(MAKE) -C ui clean

ui:
	$(MAKE) -C $@

tags:
	ctags -R
