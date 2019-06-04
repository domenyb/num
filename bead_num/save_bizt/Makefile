CFLAGS = -Wall -std=gnu99
LDLIBS = `gsl-config --libs`
#linear: plot_linear.svg lin_int.svg

#plot_linear.svg: plot.gpi pre_data.dat interpolated_data.dat Makefile
#	gnuplot $< > $@

#lin_int.svg: plot_inter.gpi interpolated_data.dat Makefile
#	gnuplot $< > $@
all: out.dat

out.dat: main
	./$<  > $@
	cat $@

main: main.c bases.h minim.h

.PHONEY: clean
clean:
	$(RM) *.o out* *.svg main out*.dat fit.log *.gch
