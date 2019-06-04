CFLAGS = -w -std=gnu99 `gsl-config --cflags` -O3
LDLIBS = $$(gsl-config --libs)

all: output
	cat $<

output: main
	./$< > $@
main: main.c extra.h

clean:
	$(RM) main output
