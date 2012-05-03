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

all: $(OTHER_PLATFORM) text

text: text.o $(DEPENDENCY_OBJ)
	$(LINKER) $^ -lm -o$@

text.o: text.c
	$(CC) -c $(CFLAGS) $^ -o$@

clean:
	-rm *.o
	-rm text win32
	-rm plot.txt
	-rm *.res *.exe *.so

wine: win32.exe

mingw: win32.exe strip

win32.exe: resource.o win32.o $(DEPENDENCY_OBJ)
	$(LINKER) -mwindows $(CFLAGS) -o$@ $^ -lgdi32 -lcomctl32

win32.o: win32.c
	$(CC) -c $(CFLAGS) $^ -o$@

strip: win32.exe
	$(STRIP) -s $^

resource.o: win32.rc
	$(RC) -o$@ $^
