.PHONY: tags all ui clean test


ifeq "$(BUILD_WITH)" "wine"
export CC:=winegcc
export RC:=wrc
endif


CFLAGS:=-O3 -Wall $(CFLAGS)

all: global.o mom.o mesh.o md.o ui

test: md_test mom_test

global.o: global.c
	$(CC) -c $(CFLAGS) $^ -o $@

mom.o: mom.c
	$(CC) -c $(CFLAGS) $^ -o $@

mesh.o: mesh.c
	$(CC) -c $(CFLAGS) $^ -o $@

md.o: md.c
	$(CC) -c $(CFLAGS) $^ -o $@

md_test: md_test.c md.c global.c
	$(CC) -D MD_ENABLE_DEBUG $(CFLAGS) $^ -lm -o $@

mom_test: mom_test.c mom.c mesh.c md.c global.c
	$(CC) -D MD_ENABLE_DEBUG -D MOM_MESH_ENABLE_DEBUG $(CFLAGS) $^ -lm -o $@

clean:
	-rm tags
	-rm *.o
	-rm mom_test md_test *.exe *.so
	-rm *.plot
	-$(MAKE) -C ui clean

ui:
	$(MAKE) -C $@

tags:
	ctags -R
