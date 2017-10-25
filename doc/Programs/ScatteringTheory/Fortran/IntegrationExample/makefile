XLF= gfortran -O3 -static

IntegrationExample :   f90lib.o IntegrationExample.o
	       ${XLF} -o  IntegrationExample  f90lib.o IntegrationExample.o 

f90lib.o :   f90lib.f90
		${XLF} -c  f90lib.f90


IntegrationExample.o :   IntegrationExample.f90
		${XLF} -c IntegrationExample.f90

