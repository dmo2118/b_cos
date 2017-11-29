CFLAGS=-Wall -O3 -g -DNDEBUG
LDLIBS=-lpng -lm

b_cos: main.o b_cos.o

clean:
	-rm *.o b_cos
