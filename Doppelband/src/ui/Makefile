.PHONY: all clean strip

LINKER:=$(CC)

CFLAGS:=-O3 -Wall $(CFLAGS)
DEPENDENCY_OBJ:= ../global.o ../mom.o ../mesh.o ../md.o
DEPENDENCY_SRC:= ../global.c ../mom.c ../mesh.c ../md.c


ifeq "$(BUILD_WITH)" "wine"
OTHER_PLATFORM:=$(OTHER_PLATFORM) wine
endif

ifeq "$(BUILD_WITH)" "mingw"
OTHER_PLATFORM:=$(OTHER_PLATFORM) mingw
endif

all: $(OTHER_PLATFORM) text$(POSTFIX)

text$(POSTFIX): text.o $(DEPENDENCY_OBJ)
	$(LINKER) $(CFLAGS) $^ -lm -o$@

text.o: text.c
	$(CC) -c $(CFLAGS) $^ -o$@

clean:
	-rm *.o
	-rm text win32
	-rm plot.txt
	-rm *.res *.exe *.so

wine: win32$(POSTFIX)
mingw: win32$(POSTFIX) strip

win32$(POSTFIX): resource.o win32_manifest.o win32.o $(DEPENDENCY_OBJ)
	$(LINKER) -mwindows $(CFLAGS) -o$@ $^ -lgdi32 -lcomctl32

win32.o: win32.c
	$(CC) -c $(CFLAGS) $^ -o$@

strip: win32$(POSTFIX) text$(POSTFIX)
	-$(STRIP) -s win32$(POSTFIX)
	-$(STRIP) -s text$(POSTFIX)

resource.o: win32.rc
	$(RC) -o$@ $^

win32_manifest.o: win32_manifest.rc
	$(RC) -o$@ $^
