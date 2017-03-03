f95 -fcheck=all -c opt_extr.f90 
f95 -fcheck=all -c marq.f90
f95 -o optimal optimal.f90 marq.o opt_extr.o -L/home/milan/cfitsio/lib -lcfitsio
rm marq.o opt_extr.o
