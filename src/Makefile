FC = gfortran -cpp -Wall -Wextra
#FC = gfortran -cpp -Wall -Wextra -g -fbounds-check

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
