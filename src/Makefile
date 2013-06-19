FC = gfortran -cpp -Wall -Wextra
#FC = gfortran -cpp -Wall -Wextra -g -fbounds-check

ho: solve.o ho.o
	$(FC) -o $@ $+

pt: solve.o pt.o
	$(FC) -o $@ $+

findpt: solve.o findpt.o
	$(FC) -o $@ $+

.PHONY: tags
tags:
	etags *.f90

test: findpt
	cd ../test && ./test.sh

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) -c $< -o $@
