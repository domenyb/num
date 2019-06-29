CFLAGS = -Wall -std=gnu99
LDLIBS = `gsl-config --libs`
#linear: plot_linear.svg lin_int.svg

#plot_linear.svg: plot.gpi pre_data.dat interpolated_data.dat Makefile
#	gnuplot $< > $@

#lin_int.svg: plot_inter.gpi interpolated_data.dat Makefile
#	gnuplot $< > $@
all: out.dat out.pdf

out.dat: main Makefile
	./$< 4 > $@
	cat $@

out.pdf: out.tex Makefile
	pdflatex $< >/dev/null

main: main.c bases.h sweep.h

.PHONEY: clean
clean:
	$(RM) *.o out.dat *.svg main times_output*.dat fit.log out.pdf out.log out.aux
