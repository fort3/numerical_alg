c= gfortran -fimplicit-none -fbounds-check -fbacktrace -g -g3 -fdefault-real-8 -O0 -finit-real=nan 

all_objs := types.o newcot.o

all: main.o $(all_objs)
	$c -g -o mynewcot main.o $(all_objs)

types.o: types.f90
	$c -c types.f90

newcot.o: types.o newcot.f90
	$c -c newcot.f90

main.o: $(all_objs) main.f90
	$c -c main.f90 

#clean:
#	rm -rf*.o *.mod
