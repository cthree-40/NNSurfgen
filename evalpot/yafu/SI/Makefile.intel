#Macro definition
FC = ifort
FFLAGS =  -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
           -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_core.a \
           $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
           -Wl,--end-group -lpthread -lm
OBJ = nn.o pes.o main.o
#end of Macro definition

test.x: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o test.x $(LDFLAGS)

clean:
	rm -f *.o *.mod a.out *.x *.exe

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
