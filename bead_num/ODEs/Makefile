CFLAGS = -Wall -std=gnu99
LDLIBS = `gsl-config --libs`
#linear: plot_linear.svg lin_int.svg

#plot_linear.svg: plot.gpi pre_data.dat interpolated_data.dat Makefile
#	gnuplot $< > $@

#lin_int.svg: plot_inter.gpi interpolated_data.dat Makefile
#	gnuplot $< > $@
all: out.dat plot.svg plot2.svg

plot.svg: plot.gpi data.dat
	gnuplot plot.gpi > plot.svg


plot2.svg: plot2.gpi data2.dat
		gnuplot plot2.gpi > plot2.svg

out.dat: main
	./$<  > $@
	cat $@

main: main.c bases.h ode.h

.PHONEY: clean
clean:
	$(RM) *.o out* *.svg main data*.dat fit.log *.gch
