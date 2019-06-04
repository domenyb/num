CFLAGS = -Wall -std=gnu99

linear: plot_linear.svg lin_int.svg

plot_linear.svg: plot.gpi pre_data.dat interpolated_data.dat Makefile
	gnuplot $< > $@

lin_int.svg: plot_inter.gpi interpolated_data.dat Makefile
	gnuplot $< > $@

pre_data.dat interpolated_data.dat: main
	./$< 2> pre_data.dat 1> interpolated_data.dat

main: main.c

.PHONEY: clean
clean:
	$(RM) *.o out* *.svg main *_data.dat
