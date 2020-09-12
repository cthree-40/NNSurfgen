#Macro definition
FC = ftn
FFLAGS = -O3
LDFLAGS =
OBJ = nn.o pes.o main.o
#end of Macro definition

test.x: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o test.x $(LDFLAGS)

clean:
	rm -f *.o *.mod a.out *.x *.exe

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
