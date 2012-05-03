.PHONY: tags all ui clean test

CFLAGS:=$(CFLAGS) -O3 -Wall

ifeq "$(BUILD_WITH)" "wine"
export CFLAGS:=$(CFLAGS)\
		-DCRT_NO_FMAX\
		-mno-cygwin\
		-DUNICODE -D_UNICODE\
		-DMOM_ENABLE_COMCTL
export CC:=winegcc
export RC:=wrc
export POSTFIX:=.exe
endif

ifeq "$(BUILD_WITH)" "mingw"
export CFLAGS:=$(CFLAGS)\
		-Dswprintf_s=snwprintf\
		-DUNICODE -D_UNICODE\
		-DMOM_ENABLE_COMCTL
export CC:=i686-w64-mingw32-gcc
export RC:=i686-w64-mingw32-windres
export STRIP:=i686-w64-mingw32-strip
export POSTFIX:=.exe
endif

all: global.o mom.o mesh.o md.o ui

test: md_test$(POSTFIX) mom_test$(POSTFIX)

global.o: global.c
	$(CC) -c $(CFLAGS) $^ -o$@

mom.o: mom.c
	$(CC) -c $(CFLAGS) $^ -o$@

mesh.o: mesh.c
	$(CC) -c $(CFLAGS) $^ -o$@

md.o: md.c
	$(CC) -c $(CFLAGS) $^ -o$@

md_test$(POSTFIX): md_test.c md.c global.c
	$(CC) -D MD_ENABLE_DEBUG $(CFLAGS) $^ -lm -o$@

mom_test$(POSTFIX): mom_test.c mom.c mesh.c md.c global.c
	$(CC) -D MD_ENABLE_DEBUG -D MOM_MESH_ENABLE_DEBUG $(CFLAGS) $^ -lm -o$@

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
