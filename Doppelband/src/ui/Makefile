.PHONY: all clean

LINKER=$(CC)

CFLAGS:=-O3 -Wall $(CFLAGS)
DEPENDENCY:= ../global.o ../mom.o ../mom_mesh.o ../md.o
DEPENDENCY_SRC:= ../global.c ../mom.c ../mom_mesh.c ../md.c

all: text win32.exe

text: text.o $(DEPENDENCY)
	$(LINKER) -lm $^ -o $@

text.o: text.c
	$(CC) -c $(CFLAGS) $^ -o $@

clean:
	-rm *.o
	-rm text
	-rm plot.txt
	-rm *.res *.exe *.exe.so
	-rm afxres.h

win32.exe: win32.c $(DEPENDENCY_SRC)
	-ln -s fake_afxres.h afxres.h
	wrc win32.rc
	winegcc $(CFLAGS) -DUNICODE -mwindows -mno-cygwin win32.res $^ -o $@
