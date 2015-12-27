all:
	mpicc -o mpi4 4mpi.c -lm
run:
	mpirun -np 2 ./mpi4 2100 59
