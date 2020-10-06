#Macro definition
FC = ifort -traceback -i8
AR = ar
FFLAGS =  -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS = $(BLAS_LIB)
OBJ = nn.o progdata.o hddata.o combinatorial.o \
		libutil.o pip.o libinternal.o localcoord.o io.o makesurf.o lm0.o lm.o
LIBO= nn.o progdata.o hddata.o combinatorial.o libutil.o pip.o libinternal.o \
		localcoord.o io.o makesurf.o lm.o potlib.o

#end of Macro definition

surfgen-nn0.x: $(OBJ) surfgen0.o main.o
	$(FC) $(FFLAGS) $(OBJ) surfgen0.o main.o -o surfgen-nn0.x $(LDFLAGS) $(CUDAFLAGS)
surfgen-nn.x: $(OBJ) surfgen.o main.o
	$(FC) $(FFLAGS) $(OBJ) surfgen.o main.o -o surfgen-nn.x $(LDFLAGS) $(CUDAFLAGS)
	cp surfgen-nn.x ../bin/oh3_surfgen-nn.cinput.x
libs: $(LIBO)
	$(AR) -r -v lib-nnsurfgen.a $(LIBO)
	cp lib-nnsurfgen.a ../lib/
evalpot: $(LIBO) evaluatepot.o
	$(FC) $(FFLAGS) evaluatepot.o -o evalpot.x $(LDFLAGS) ../lib/lib-nnsurfgen.a
clean:
	rm -f *.o *.mod a.out *.x *.a *.exe

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
