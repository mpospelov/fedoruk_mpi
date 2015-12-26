#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include <mpi.h>

int main(int argc, char **argv) {
  int myrank, total;

  double *A, *B, *C;	// Используются только в root

  double *a, *b, *c;	// Лента матрицы [mxn], вектор [n], рез-т [m]

  int n;    // Размерность квадратной матрицы
  int m;    // Ширина горизонтальной ленты матрицы
  int i, j;
  int intBuf[2];

  MPI_Init (&argc, &argv); // Инициализация коммуникационных средств
  MPI_Comm_size (MPI_COMM_WORLD, &total);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  printf ("Total=%d, rank=%d\n", total, myrank);

  if (!myrank) {
    // Подготовка исх. данных (только root)
    n = N_PER_PROC * total;
    A = (double *) malloc (sizeof(double)*n*n);
    B = (double *) malloc (sizeof(double)*n);
    C = (double *) malloc (sizeof(double)*n);
    // Инициализация матрицы A и вектора B
    for (i=0; i<n; i++) {
      B[i] = (double)i;
      for (j=0; j<n; j++)
        A[i*n+j] = (double)(i+j);
    }
  }

  if (!myrank) {
    intBuf[0] = n;
    intBuf[1] = N_PER_PROC;
  }
  MPI_Bcast((void *)intBuf, 2, MPI_INT, 0, MPI_COMM_WORLD); // Широковещательная рассылка/прием сообщения с блокировкой процессов

  n = intBuf[0];
  m = intBuf[1];

  a = (double *) malloc (sizeof(double)*n*m);
  b = (double *) malloc (sizeof(double)*n);
  c = (double *) malloc (sizeof(double)*m);

  if (!myrank) {	// Лишнее действие, B не нужен
    memcpy (b, B, sizeof(double)*n);
  }
  MPI_Bcast((void *)b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Широковещательная рассылка/прием

  MPI_Scatter((void *)A, n*m, MPI_DOUBLE,
	      (void *)a, n*m, MPI_DOUBLE, 0, MPI_COMM_WORLD); // одновременная рассылка разных (но однотипных) данных разным процессам в группе

  for (i=0; i<m; i++) {
    c[i] = 0;
    for (j=0; j<n; j++)
      c[i] += a[n*i+j]*b[j];
  }

  MPI_Gather((void *)c, m, MPI_DOUBLE,
	    (void *)C, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

 /* if (!myrank)
    for (i=0; i<n; i++)
      printf ("%g\n", C[i]); */

  MPI_Finalize(); // Нормальное завершение обменов

  exit(0);
}
