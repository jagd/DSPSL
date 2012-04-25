.PHONY: all clean

LINKER=$(CC)

CFLAGS:=-O3 -Wall $(CFLAGS)
DEPENDENCY:= ../global.o ../mom.o ../mesh.o ../md.o
DEPENDENCY_SRC:= ../global.c ../mom.c ../mesh.c ../md.c

all: text win32.exe
	$(MAKE) text
	-$(MAKE) win32.exe

text: text.o $(DEPENDENCY)
	$(LINKER) -lm $^ -o $@

text.o: text.c
	$(CC) -c $(CFLAGS) $^ -o $@

clean:
	-rm *.o
	-rm text
	-rm plot.txt
	-rm *.res *.exe *.exe.so

win32.exe: win32.rc win32.c $(DEPENDENCY_SRC)
	wrc win32.rc
	winegcc $(CFLAGS) -DUNICODE -mwindows -mno-cygwin \
		win32.res win32.c\
		$(DEPENDENCY_SRC) -o $@
