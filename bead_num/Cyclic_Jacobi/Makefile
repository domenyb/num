CFLAGS = -Wall -std=gnu99
LDLIBS = `gsl-config --libs`
#linear: plot_linear.svg lin_int.svg

#plot_linear.svg: plot.gpi pre_data.dat interpolated_data.dat Makefile
#	gnuplot $< > $@

#lin_int.svg: plot_inter.gpi interpolated_data.dat Makefile
#	gnuplot $< > $@
all: out.dat out.svg

out.svg: plt.gpi times_output.dat
		gnuplot $< > $@ 2> /dev/null

times_output.dat: main Makefile
	for i in `seq 100 2 200`; do \time -f "$$i %U" -ao $@ ./main $$i > /dev/null ; done

#	for i in `seq 5 2 121`; do time --format "$$i %U" -ao time.txt $@ > /dev/null; done

out.dat: main
	./$< 5 > $@
	cat $@

main: main.c bases.h sweep.h

.PHONEY: clean
clean:
	$(RM) *.o out* *.svg main times_output.dat fit.log
