
FC=gfortran

FFLAGS=-w -fno-automatic -ffixed-line-length-132
#FFLAGS=-fPIC -w -fno-automatic -ffixed-line-length-132

obj=Func_PSI.o genps_nbody.o

all: libgenps.a

libgenps.a: $(obj)
	@rm -f libgenps.a
	@ar cr libgenps.a $(obj)
	@rm *.o

test: test.o $(obj)
	$(FC) $(FFLAGS) -o $@ test.o $(obj)

test.o: test.f90 genps_nbody.mod
	$(FC) $(FFLAGS) -c $< -o $@ -lstdc++
#	$(FC) -c $< $(FFLAGS)

Func_PSI.o Func_PSI.mod: Func_PSI.f90
	$(FC) $(FFLAGS) -c $< -o $@ -lstdc++
#	$(FC) -c $< $(FFLAGS)


genps_nbody.o genps_nbody.mod: genps_nbody.f90 Func_PSI.mod
	$(FC) $(FFLAGS) -c $< -o $@ -lstdc++
#	$(FC) -c $< $(FFLAGS)

clean:
	rm *~ *.o *.mod test libgenps.a
