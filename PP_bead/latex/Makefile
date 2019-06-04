CFLAGS = -w -std=gnu99 `gsl-config --cflags`
LDLIBS = $$(gsl-config --libs)

default: out.dat error.pdf Makefile
	cat $<

error.pdf: error.tex error_plot.png
	pdflatex $<


error_plot.png: plot.gpi plot.data
	gnuplot plot.gpi > error_plot.png

plot.data:main
	./main -3 3 0.2  > out.dat

clean:
	$(RM) *.o out* *.svg main plot.data *.pdf *.log *.aux plot-cairo.tex
