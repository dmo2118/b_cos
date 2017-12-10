# GPROF=-pg

CFLAGS=-Wall -g $(GPROF) -O3 -DNDEBUG
LDLIBS=-lpng -lm $(GPROF)

b_cos: main.o b_cos.o

clean:
	-rm *.o b_cos
