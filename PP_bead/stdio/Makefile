CFLAGS = -Wall -std=gnu99
LDLIBS = -lm


all:out.io.txt out.cmd.dat plot.svg

out.io.txt:main-stdio input_data
	echo 1 2 3 4 5 | ./main-stdio > out.io.txt
	cat out.io.txt
	seq 1 0.5 4 | sed -e 's/,/./g' | ./main-stdio > io.out.txt
	cat out.io.txt
	cat input_data | ./main-stdio >io.out.txt

out.cmd.dat: main-cmdline input_data
	./main-cmdline 1 2 3 4 5 6 > out.cmd.dat
	cat out.cmd.dat
	./main-cmdline `seq 0 0.5 2` | sed -e 's/,/./g' > out.cmd.dat
	cat out.cmd.dat
	./main-cmdline `cat input_data` > out.cmd.dat

input_data:
	seq 0 0.2 5 | sed -e 's/,/./g' > input_data

plot.svg:plot.gpi
	gnuplot plot.gpi

main-stdio: main-stdio.c

main-cmdline: main-cmdline.c

clean:
	$(RM) main-stdio main-cmdline *.txt *.dat plot.svg input_data
