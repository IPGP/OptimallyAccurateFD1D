EXE= FWI
OBJS=Opt2.o main.o data.o pml.o
FC = gfortran								
FFLAGS= -O3 -Wall -o


all:Opt2 data  pml main FWI

#main file from where all the functions are launched	
main:main.f95    
	$(FC) -Wall -c main.f95 

# Forward modelling with optimally accurate operators
Opt2:Opt2.f95

	$(FC) -Wall -c Opt2.f95

# data
data:data.f95
	$(FC)  -Wall -c data.f95 

# pml sigma
pml:pml.f95
	$(FC)  -Wall -c  pml.f95

#exectutable file (all the above files are linked)
FWI:Opt2.o main.o data.o pml.o

	$(FC)  -Wall -o FWI Opt2.o pml.o main.o data.o 


clean:
	rm $(EXE) $(OBJS)	
