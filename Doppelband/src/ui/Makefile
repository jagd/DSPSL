.PHONY: all clean

LINKER=$(CC)

CFLAGS:=-O3 -Wall $(CFLAGS)
DEPENDENCY:= ../mom.o ../mom_mesh.o ../md.o

all: text

text: text.o $(DEPENDENCY)
	$(CC) -lm $^ -o $@

text.o: text.c
	$(LINKER) -c $(CFLAGS) $^ -o $@

clean:
	-rm *.o
	-rm text
	-rm plot.txt
