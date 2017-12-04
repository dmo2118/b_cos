# GPROF=-pg

CFLAGS=-Wall -g -DNDEBUG $(GPROF) # -O3
LDLIBS=-lpng -lm $(GPROF)

b_cos: main.o b_cos.o

clean:
	-rm *.o b_cos
