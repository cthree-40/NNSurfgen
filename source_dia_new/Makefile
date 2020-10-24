#Macro definition
FC = ifort -traceback -i8
FFLAGS = -g -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS = $(BLAS_LIB)
OBJ = nn.o progdata.o hddata.o pip.o libinternal.o \
	io.o makesurf.o lm.o surfgen.o main.o
#end of Macro definition

surfgen-nn.x: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o surfgen-nn-dia.x $(LDFLAGS) $(CUDAFLAGS)
	cp surfgen-nn-dia.x ../bin/oh3_surfgen-nn-dia.x
clean:
	rm -f *.o *.mod a.out *.x *.a *.exe

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
