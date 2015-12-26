all:
	mpicc -o mpi4 4mpi.c -lm
run:
	mpirun -np 3 ./mpi4 1500 8000
