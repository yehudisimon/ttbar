
FC=gfortran

FFLAGS=-w -fno-automatic -ffixed-line-length-132
#FFLAGS=-fPIC -w -fno-automatic -ffixed-line-length-132

objects=nintlib.o

all: libnin.a

nintlib.o nintlib.mod: nintlib.f90
	$(FC) $(FFLAGS) -c $< -o $@ -lstdc++

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@ -lstdc++

libnin.a: $(objects)
	@rm -f libnin.a
	@ar cr libnin.a $(objects)
	@rm *.o

clean:
	@rm libnin.a *~ *.o
