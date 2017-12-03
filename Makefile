# GPROF=-pg

CFLAGS=-Wall -O3 -g -DNDEBUG $(GPROF)
LDLIBS=-lpng -lm $(GPROF)

b_cos: main.o b_cos.o

clean:
	-rm *.o b_cos
