CFLAGS = -pthread -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)


all: out.txt Makefile
		cat out.txt

out.txt plot.dat:main
		./$< 1> out.txt 2> plot.dat

clean:
	$(RM) main main.o out.txt plot.svg plot.dat
