CFLAGS = -w -std=gnu99 `gsl-config --cflags`
LDLIBS = $$(gsl-config --libs)

all: output
	cat $<

output: main
	./$< > $@

clean:
	$(RM) main output
