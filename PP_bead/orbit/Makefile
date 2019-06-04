CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)


all:plot.svg plot2.svg


plot.svg:plot.gpi plot.dat
	gnuplot $< > $@

plot2.svg:plot2.gpi plot2.dat
		gnuplot $< > $@

plot.dat plot2.dat:main
	./main 1> plot.dat 2> plot2.dat


clean:
	$(RM) *.o out* *.svg main plot*.dat
