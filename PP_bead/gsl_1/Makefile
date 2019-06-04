CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = `gsl-config --libs`

all: airy.svg out.txt
	cat out.txt

airy.svg: airy.gpi airy.dat
	gnuplot airy.gpi > $@

out.txt airy.dat:main
	cat input_data | ./$< 1> out.txt 2> airy.dat

main: main.c airy.h


clean:
	$(RM) *.o out* *.svg main airy.dat
