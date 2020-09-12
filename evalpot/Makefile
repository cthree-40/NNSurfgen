#Macro definition
FC = ftn -traceback -i8
FFLAGS = -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS = $(BLAS_LIB)
EVALPOTOBJ = nn.o pes.o evalpot.o
TESTGRADOBJ = nn.o pes.o testgrads.o
TESTPTSOBJ = nn.o pes.o testpoints.o
FINDCPOBJ  = nn.o pes.o opttools.o findcp.o
FINDMXOBJ  = nn.o pes.o opttools.o findmex.o
#end of Macro definition

evalpot.x: $(EVALPOTOBJ)
	$(FC) $(FFLAGS) $(EVALPOTOBJ) -o evalpot.x $(LDFLAGS)

findcp.x:  $(FINDCPOBJ)
	$(FC) $(FFLAGS) $(FINDCPOBJ) -o findcp.x $(LDFLAGS)

findmex.x:  $(FINDMXOBJ)
	$(FC) $(FFLAGS) $(FINDMXOBJ) -o findmex.x $(LDFLAGS)

testgrads.x: $(TESTGRADOBJ)
	$(FC) $(FFLAGS) $(TESTGRADOBJ) -o testgrads.x $(LDFLAGS)

testpoints.x: $(TESTPTSOBJ)
	$(FC) $(FFLAGS) $(TESTPTSOBJ) -o testpoints.x $(LDFLAGS)

clean:
	rm -f *.o *.mod a.out *.x *.exe

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
