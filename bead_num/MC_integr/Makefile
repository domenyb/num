CFLAGS = -Wall -std=gnu99
LDLIBS = `gsl-config --libs`
#linear: plot_linear.svg lin_int.svg

#plot_linear.svg: plot.gpi pre_data.dat interpolated_data.dat Makefile
#	gnuplot $< > $@

#lin_int.svg: plot_inter.gpi interpolated_data.dat Makefile
#	gnuplot $< > $@
all: out.dat out.svg

out.svg: plot.gpi data.dat
	gnuplot $< > $@

out.dat data.dat: main
	./$<  > out.dat
	cat out.dat

main: main.c bases.h integration.h

.PHONEY: clean
clean:
	$(RM) *.o out* *.svg main data*.dat fit.log *.gch
