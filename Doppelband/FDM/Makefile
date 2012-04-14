.PHONY: all clean

CFLAGS:=-O3 -Wall $(CFLAGS)

all: CSR_Test

CSR_Test: CSR_Test.c CSR.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	-rm *.o
	-rm CSR_Test
