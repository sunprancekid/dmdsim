#!/user/bin/env make

# specify the fortran compiler
FORT:=gfortran

# location of simclass files
BIN=./fortran

# debugging / testing flags
DEBUG=-fcheck=all -fbounds-check


default:
	echo "No command given."

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exc
	rm -f *.dat
	rm -f *.xyz
	rm -f a.out
	rm -f fort.*
	
open: 
	sublime *.f90
	sublime ${BIN}/*.f90


# active particle program
CONSTANT=$(BIN)/dmdconstant
BINARYTREE=$(BIN)/binarytree
TYPEID=$(BIN)/type_id
ACTIVEDISK=actdisk

constant: 
	$(FORT) -c $(CONSTANT).f90 -o $(CONSTANT).o

binarytree: 
	$(FORT) -c $(BINARYTREE).f90 -o $(BINARYTREE).o

type_id: constant
	$(FORT) -c $(TYPEID).f90 -o $(TYPEID).o

act: constant type_id binarytree
	$(FORT) -c $(ACTIVEDISK).f90 -o $(ACTIVEDISK).o
	
act_test: clean  act
	$(FORT) $(DEBUG) -O3 $(CONSTANT).o $(BINARYTREE).o $(TYPEID).o $(ACTIVEDISK).o
	./a.out
	


