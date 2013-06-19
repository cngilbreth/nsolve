FC = gfortran -cpp -Wall -Wextra
#FC = gfortran -cpp -Wall -Wextra -g -fbounds-check

vpath ho ../bin
vpath pt ../bin

ho: solve.o ho.o
	$(FC) -o ../bin/$@ $+

pt: solve.o pt.o
	$(FC) -o ../bin/$@ $+

findpt: solve.o findpt.o
	$(FC) -o ../bin/$@ $+

.PHONY: check
check: ho pt
	@cd ../test && ./testho.sh && ./testpt.sh

.PHONY: tags
tags:
	etags *.f90

clean:
	rm -f *.o *.mod

%.o: %.f90
	$(FC) -c $< -o $@
