CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)


all: out.txt Makefile plot.svg
		cat out.txt

plot.svg: plot.gpi plot.dat
		gnuplot $< > $@

out.txt plot.dat:main
		./$< 1> out.txt 2> plot.dat

clean:
	$(RM) main main.o out.txt plot.svg plot.dat
