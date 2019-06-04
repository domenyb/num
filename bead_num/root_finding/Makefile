CFLAGS = -Wall -std=gnu99
LDLIBS = `gsl-config --libs`
#linear: plot_linear.svg lin_int.svg

#plot_linear.svg: plot.gpi pre_data.dat interpolated_data.dat Makefile
#	gnuplot $< > $@

#lin_int.svg: plot_inter.gpi interpolated_data.dat Makefile
#	gnuplot $< > $@
all: out.dat out.svg

out.svg: plt.gpi
		gnuplot $< > $@

out.dat: main input_vectors.dat
	./$< input_vectors.dat > $@
	cat $@

main: main.c bases.h root.h

.PHONEY: clean
clean:
	$(RM) *.o out* *.svg main out*.dat fit.log *.gch
