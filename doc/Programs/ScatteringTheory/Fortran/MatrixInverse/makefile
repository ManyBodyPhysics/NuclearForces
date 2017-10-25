XLF= gfortran -O3 -static

MatrixInverse :   f90lib.o MatrixInverse.o
	       ${XLF} -o  MatrixInverse  f90lib.o MatrixInverse.o 

f90lib.o :   f90lib.f90
		${XLF} -c  f90lib.f90


MatrixInverse.o :   MatrixInverse.f90
		${XLF} -c MatrixInverse.f90

