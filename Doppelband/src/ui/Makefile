.PHONY: all clean

LINKER=$(CC)

CFLAGS:=-O3 -Wall $(CFLAGS)
DEPENDENCY_OBJ:= ../global.o ../mom.o ../mesh.o ../md.o
DEPENDENCY_SRC:= ../global.c ../mom.c ../mesh.c ../md.c


ifeq "$(BUILD_WITH)" "wine"
OTHER_PLATFORM:=$(OTHER_PLATFORM) wine
endif

ifeq "$(BUILD_WITH)" "mingw"
OTHER_PLATFORM:=$(OTHER_PLATFORM) mingw
endif

all: text $(OTHER_PLATFORM)

text: text.o $(DEPENDENCY_OBJ)
	$(LINKER) $^ -lm -o $@

text.o: text.c
	$(CC) -c $(CFLAGS) $^ -o $@

clean:
	-rm *.o
	-rm text win32
	-rm plot.txt
	-rm *.res *.exe *.so

wine: win32

mingw: win32

win32: win32.res win32.o $(DEPENDENCY_OBJ)
	$(LINKER) -DUNICODE -mwindows -mno-cygwin -lcomctl32 $^ -o $@

win32.o: win32.c
	$(CC) -c -D MOM_ENABLE_COMCTL $(CFLAGS) $^ -o $@

win32.res: win32.rc
	$(RC) win32.rc
