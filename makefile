CXX=g++-11
CC=gcc-11
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas -lfftw3 -lblas

CDEPS = AYaux.h

CPPDEPS = AYlinalg.hh

OBJS1=main1.o AYaux.o AYlinalg.o AYmat.o AYvec.o AYrng.o AYtens.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

all: test1

test1: $(OBJS1)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

clean:
	rm -f test1
	rm -f *.o
clean_dat:
	rm -f *.dat
	rm -f *.aysml
	rm -f *.aydat
clean_datdir:
	rm -f ./dat_dir/*.dat
	rm -f ./dat_dir/*.aysml
	rm -f ./dat_dir/*.aydat

clean_all: clean clean_dat
