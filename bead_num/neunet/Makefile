CFLAGS = -Wall -std=gnu99
LDLIBS = `gsl-config --libs`
#linear: plot_linear.svg lin_int.svg

#plot_linear.svg: plot.gpi pre_data.dat interpolated_data.dat Makefile
#	gnuplot $< > $@

#lin_int.svg: plot_inter.gpi interpolated_data.dat Makefile
#	gnuplot $< > $@
all: out.dat plot.svg plot2.svg plot3.svg

plot2.svg: plot2.gpi data2.dat
		gnuplot $< > $@

plot3.svg: plot3.gpi data3.dat
			gnuplot $< > $@


plot.svg: plot.gpi data.dat
	gnuplot $< > $@

out.dat data.dat data2.dat data3.dat: main
	./$<  > out.dat
	cat out.dat

main: main.c ann.h

.PHONEY: clean
clean:
	$(RM) *.o out* *.svg main data*.dat fit.log *.gch
