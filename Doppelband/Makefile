.PHONY: all clean

CFLAGS:=-O3 -Wall $(CFLAGS)

all: md_test

md_test: md_test.c md.c
	$(CC) -D MD_ENABLE_DEBUG  $(CFLAGS) $^ -o $@

clean:
	-rm *.o
	-rm md_test
