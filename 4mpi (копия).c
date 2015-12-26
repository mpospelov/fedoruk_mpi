#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include <mpi.h>
//
#include <stdio.h>  
#include <string.h>  
#include <math.h>

// mpicc -o mpi4 mpi4.c -lm
// mpirun ./mpi4

#define debug 1 // 1 - on, 0 - off
#define gnu 0 // 1 - on, 0 - off
#define flagForPrintArrayTOGNUANDSOLV 0

#define defCondOne -999999 // inicialize first border condition

#define _REENTRANT

#define N_PER_PROC 100  // Кол-во строк матрицы на процесс

double calcTime;
double dx, dy, dt;
double length, hight;
long xPoints, yPoints;
long yPointsMPI;
long sendRowsMPI;
double at;

double printDx, printDy;
long printXPoints, printYPoints;

long numPoints;

double maxT;
double minT;

FILE *file;

double** tempMas;
double** tempMasPrev;

double** tempMas2;
double** tempMasPrev2;

double** borderConditionOne;
double** borderConditionTwo;
 
double** make2DArray(long arraySizeX, long arraySizeY);

void inicialBeginTemp(double **, double);

void printMas(double**, long, long);
void printMasI(double*, long, long);
long R(long, long, long, long );
void solve(double** , double** , double** , double** , long);
void inicialBorderConditionsOne( double** , double );
void inicialBorderConditionsTwo( double** , double );
void applyBorderOneToTemp(double** ,double**);

void makeArrayForSendToSolve(double*, double**, int);
void makeArrayForSendToGNU(double*, double**, int);
void MPIsolveBorderOneTwo(double** , double**, double** , double**); 
void solveMPI(double* , double*, double** , double**, long, long);

void findMinMaxGnuPlot(double** );
void makeFinalFileFileOpen();
void makeFinalFileSETUP_phase0(double** );
void makeFinalFileDATA_phase1(double** );
void makeFinalFileEND_phase2();

void compareTwoArray(double**, double**, long, long);

int main( int argc, char **argv ){

  dt = 0.1;

  length = 2.0;
  hight = 1.1;

  if (argc < 3)
  {
    printf("\nError input parameters!\n");
      
    calcTime = 1500.01; // seconds
    //calcTime = 0.01; // seconds

    dx = 0.4;
    dy = 0.1;

  } else {
    calcTime = atoi(argv[1]);
    numPoints = atoi(argv[2]);

    int x = sqrt( numPoints / ((length) * (hight)) );
    printf("\n x= %d", x);
    dx = (double)(length / ( ((length) * x) -1.0) );
    dy = (double)(hight / ( ((hight) * x) -1.0) );
    printf("\n dx= %f(%f), dy= %f(%f) ", 
      (double)dx, (double)dx*(length*((int)x) -1.0), 
      (double)dy, (double)dy*(hight*((int)x) -1.0) );

    printf("\n printf dx= %ld, print dy= %ld", (long)(length/ dx), (long)(hight/ dy));
  }
 
  at = 420.0 / (235.0 * 10500.0);

  xPoints = 1 + length/dx;
  yPoints = 1 + hight/dy;

  tempMas = make2DArray(xPoints, yPoints);
  tempMasPrev = make2DArray(xPoints, yPoints); 

  tempMas2 = make2DArray(xPoints, yPoints);
  tempMasPrev2 = make2DArray(xPoints, yPoints);

  inicialBeginTemp(tempMas, 0.0);
  inicialBeginTemp(tempMasPrev, 0.0);
  inicialBeginTemp(tempMas2, 0.0);
  inicialBeginTemp(tempMasPrev2, 0.0);
  int i2,j2;
  for (i2 = 0; i2 < yPoints; i2++)
  {
    for (j2 = 0; j2 < xPoints; j2++)
    {
      tempMas[i2][j2] = i2+1;
      tempMasPrev[i2][j2] = i2+1;
      tempMas2[i2][j2] = i2+1;
      tempMasPrev2[i2][j2] = i2+1;
      if (1) {
        tempMas[i2][j2] = 0.0;
        tempMasPrev[i2][j2] = 0.0;
        tempMas2[i2][j2] = 0.0;
        tempMasPrev2[i2][j2] = 0.0;
      }
    }
  }


  borderConditionOne = make2DArray(xPoints, yPoints); 
  borderConditionTwo = make2DArray(xPoints, yPoints); 

  inicialBeginTemp(borderConditionOne, defCondOne);
  inicialBeginTemp(borderConditionTwo, defCondOne);

  inicialBorderConditionsOne(borderConditionOne, 100.0); // temp in x[i][0] = temp
  inicialBorderConditionsTwo(borderConditionTwo, 300.0); // heat temp in x[i][last] = temp

  if (1) {
    printMas( tempMas, xPoints, yPoints );
    printMas( tempMasPrev, xPoints, yPoints );
    printMas( borderConditionOne, xPoints, yPoints ); 
    printMas( borderConditionTwo, xPoints, yPoints ); 
  } 

  printf("\n\n START SOLVE \n\n");

  applyBorderOneToTemp(borderConditionOne, tempMasPrev);
 
  if (gnu == 1) {
    findMinMaxGnuPlot(tempMasPrev); 
    makeFinalFileFileOpen();
  }

  // MPI
  int myrank, total;  

  double* resizeTempMas;
  double* resizeTempMasPrev;
  double* sendT;
  double* sendTPrev;

  MPI_Init (&argc, &argv); // Инициализация коммуникационных средств
  MPI_Comm_size (MPI_COMM_WORLD, &total);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  
  printf ("Total=%d, rank=%d\n", total, myrank);

  sendRowsMPI = (yPoints / total) + 2;
  yPointsMPI = yPoints + 2*total;
 
  //sendT = make2DArray(xPoints, sendRowsMPI);
  //sendTPrev = make2DArray(xPoints, sendRowsMPI);
  sendT = (double*) malloc(xPoints*sendRowsMPI*sizeof(double));
  sendTPrev = (double*) malloc(xPoints*sendRowsMPI*sizeof(double));
  for (i2 = 0; i2 < xPoints*sendRowsMPI; i2++) {
    sendT[i2] = 0.0;
    sendTPrev[i2] = 0.0;
  }

  printf("\n num= %d, begin size= %d, end size= %d \n", 
    myrank, (int)yPoints, (int)sendRowsMPI); 

  if (myrank == 0) { 
    resizeTempMas = (double*) malloc(xPoints*yPointsMPI*sizeof(double));
    resizeTempMasPrev = (double*) malloc(xPoints*yPointsMPI*sizeof(double));
  }
   
  long i;
  int flag = 0;
  for (i = 0; i < (calcTime / dt); i++) {
  
    if (flag == 0) {

      if (myrank == 0) { 
        MPIsolveBorderOneTwo(tempMas, tempMasPrev, borderConditionOne, borderConditionTwo); 

        makeArrayForSendToSolve(resizeTempMas, tempMas, total);

        solve(tempMas2, tempMasPrev, borderConditionOne, borderConditionTwo, i);
      }

      MPI_Scatter((void *)resizeTempMas, sendRowsMPI*xPoints, MPI_DOUBLE,
        (void *)sendTPrev, sendRowsMPI*xPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD); // одновременная рассылка разных (но однотипных) данных разным процессам в группе

      solveMPI(sendT, sendTPrev, borderConditionOne, borderConditionTwo, i, myrank);

      MPI_Barrier(MPI_COMM_WORLD);
      printf("\nGather MPI\n");
      MPI_Gather((void *)sendT, sendRowsMPI*xPoints, MPI_DOUBLE,
          (void *)resizeTempMasPrev, sendRowsMPI*xPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      printf("\nGather End MPI\n");

       

      if (myrank == 0) {
        printMasI(resizeTempMasPrev, xPoints, yPointsMPI);

        printf("\nMake array for GNU\n");
        makeArrayForSendToGNU(resizeTempMasPrev, tempMasPrev, total);

        compareTwoArray(tempMasPrev, tempMas2, xPoints, yPoints);
      }

      flag =1;
    } else {

      if (myrank == 0) { 
        MPIsolveBorderOneTwo(tempMasPrev, tempMas, borderConditionOne, borderConditionTwo);  

        makeArrayForSendToSolve(resizeTempMasPrev, tempMasPrev, total);
        solve(tempMasPrev2, tempMas, borderConditionOne, borderConditionTwo, i);
      }

      MPI_Scatter((void *)resizeTempMasPrev, sendRowsMPI*xPoints, MPI_DOUBLE,
        (void *)sendT, sendRowsMPI*xPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD); // одновременная рассылка разных (но однотипных) данных разным процессам в группе

      solveMPI(sendTPrev, sendT, borderConditionOne, borderConditionTwo, i, myrank);

      MPI_Barrier(MPI_COMM_WORLD);
      printf("\nGather MPI\n");
      MPI_Gather((void *)sendTPrev, sendRowsMPI*xPoints, MPI_DOUBLE,
          (void *)resizeTempMas, sendRowsMPI*xPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      printf("\nGather End MPI\n");

      if (myrank == 0) {
        printMasI(resizeTempMas, xPoints, yPointsMPI);

        printf("\nMake array for GNU\n");
        makeArrayForSendToGNU(resizeTempMas, tempMas, total);

        compareTwoArray(tempMas, tempMasPrev2, xPoints, yPoints);
      }

      flag =0;
    }
  }

 /*if (myrank == -1) {
    long i;
    int flag = 0;
    for (i = 0; i <= (calcTime / dt); i++) {

      if (myrank == 0) {
        makeArrayForSendToSolve(resizeTempMas, tempMas, total);
      }

      if (flag == 0) {
        MPI_Scatter((void *)resizeTempMas, (yPoints/total)+2, MPI_DOUBLE,
        (void *)sendT, (yPoints/total)+2, MPI_DOUBLE, 0, MPI_COMM_WORLD); // одновременная рассылка разных (но однотипных) данных разным процессам в группе

        solveMPI(sendT, sendTPrev, borderConditionOne, borderConditionTwo, i);

        flag = 1;
      } else {
        MPI_Scatter((void *)resizeTempMas, (yPoints/total)+2, MPI_DOUBLE,
        (void *)sendTPrev, (yPoints/total)+2, MPI_DOUBLE, 0, MPI_COMM_WORLD); // одновременная рассылка разных (но однотипных) данных разным процессам в группе

        solveMPI(sendTPrev, sendT, borderConditionOne, borderConditionTwo, i);

        flag = 0;
      } 
    }

    if (gnu == 1) {
      makeFinalFileEND_phase2();  
    }
  }*/


  /////////////// MPI
 
  if (myrank == -1) {
    long i;
    int flag = 0;
    for (i = 0; i <= (calcTime / dt); i++) {
      if (flag == 0) {
        solve(tempMas, tempMasPrev, borderConditionOne, borderConditionTwo, i); 
        flag = 1;
      } else {
        solve(tempMasPrev, tempMas, borderConditionOne, borderConditionTwo, i);
        flag = 0;
      } 
    }

    if (gnu == 1) {
      makeFinalFileEND_phase2();  
    }
  }

  MPI_Finalize(); // Нормальное завершение обменов

  return 0;
}

long R(long i, long k, long min, long max) {
  long ans = i + k;
  //printf("max=%ld\n", max);
  if ((ans < max) && (ans > min)) {
    //printf("ans=%ld\n",ans);
    return ans;
  
  } else {

    if (ans >= max) {
      //printf("ans=%ld\n",ans);
      return max - 1;
    } else {
      //printf("ans=%ld\n",ans);
      return min;
    }
  }
}

long IR(long i, long myrank) {
  return (i + sendRowsMPI * myrank) -1 -2*myrank;
}

void solveMPI(double* temp, double* tempPrev,
  double** borderOne, double** borderTwo, long calcIterat, 
  long myrank) {

  printf("\nSolve MPI № %ld, from %ld to %ld \n", myrank, 
    IR(1, myrank), IR(sendRowsMPI-2, myrank) );

  long i, j, k;  

  for (i = 0; i < sendRowsMPI; i++) {
    for (j = 0; j < xPoints; j++){
      printf("%6.2f ", tempPrev[i*xPoints +j]);
    }
    printf("\n ");
  }

/*
  for (i = 1; i < sendRowsMPI-1; i++) {
    for (j = 0; j < xPoints; j++){

      if (borderTwo[IR(i, myrank)][j] != defCondOne) {
        if (j == xPoints -1) {
          temp[i*xPoints +j] = borderTwo[IR(i, myrank)][j] * dx + tempPrev[i*xPoints +j-1];
          //if (i == 0) printf("%f,", temp[i][j]);
        }
      }
    }
  }*/
 
  for (i = 1; i < sendRowsMPI-1; i++) {
    for (j = 0; j < xPoints; j++){ 
      if (borderTwo[IR(i, myrank)][j] == defCondOne)
      {
        if (borderOne[IR(i, myrank)][j] == defCondOne) { // without first border condition

          double tdx, tdy;
          
          tdy = (double) ( tempPrev[i*xPoints +R(j, 1, 0, xPoints)] 
            - 2.0 * tempPrev[i*xPoints +j] 
            + tempPrev[i*xPoints +R(j, -1, 0, xPoints)] ) / (double) ( dx * dx );

          tdx = (double) ( tempPrev[R(i, 1, 0, sendRowsMPI)*xPoints +j] 
            - 2.0 * tempPrev[i*xPoints +j] 
            + tempPrev[R(i, -1, 0, sendRowsMPI)*xPoints +j] ) / (double) ( dy * dy );

          temp[i*xPoints +j] = (double) dt * at * ( tdx + tdy ) + tempPrev[i*xPoints +j];
        } else {
          temp[i*xPoints +j] = borderOne[IR(i, myrank)][j];
        }
      }
    } 
  } 
  printf("\nEnd Solve MPI № %ld \n", myrank);
}

void MPIsolveBorderOneTwo(double** temp, double** tempPrev,
  double** borderOne, double** borderTwo) {
  long i, j, k; 


  for (i = 0; i < yPoints; i++) {
    for (j = 0; j < xPoints; j++){

      if (borderTwo[i][j] != defCondOne) {
        if (j == xPoints -1) {
          temp[i][j] = borderTwo[i][j] * dx + tempPrev[i][j-1];
          //if (i == 0) printf("%f,", temp[i][j]);
        }
      } else {
        if (borderOne[i][j] != defCondOne) {
          temp[i][j] = borderOne[i][j];
        }
      } 
    }
  }
}

void solve(double** temp, double** tempPrev, 
  double** borderOne, double** borderTwo, long calcIterat) {

  long i, j, k; 

  for (i = 0; i < yPoints; i++) {
    for (j = 0; j < xPoints; j++){

      if (borderTwo[i][j] != defCondOne) {
        if (j == xPoints -1) {
          temp[i][j] = borderTwo[i][j] * dx + tempPrev[i][j-1];
          //if (i == 0) printf("%f,", temp[i][j]);
        }
      }
    }
  }

  for (i = 0; i < yPoints; i++) {
    for (j = 0; j < xPoints; j++){ 
      if (borderTwo[i][j] == defCondOne)
      {
        if (borderOne[i][j] == defCondOne) { // without first border condition

          double tdx1,tdx2, tdy1,tdy2;
          double temNull = tempPrev[i][j];

          if (j == 0) { 
            tdx1 = temNull;

          } else {
            tdx1 = tempPrev[i][ j -1];
          }

          if (j == xPoints -1) {
            tdx2 = temNull;
          } else {
            tdx2 = tempPrev[i][ j +1];
          }

          if (i == 0) {
            tdy1 = temNull;
          } else {
            tdy1 = tempPrev[ i -1][j];
          }

          if (i == yPoints -1) {
            tdy2 = temNull;
          } else {
            tdy2 = tempPrev[ i +1][j];
          }

          double tdx, tdy;
          
          tdy = (double) ( tdy2 -2.0 * tempPrev[i][j] +tdy1 ) / 
            (double) ( dy * dy ); 

          tdx = (double) ( tdx2 -2.0 * tempPrev[i][j] +tdx1 ) / 
            (double) ( dx * dx );

          temp[i][j] = (double) dt * at * ( tdx + tdy ) + tempPrev[i][j];

        } else { // first border condition
          
          temp[i][j] = borderOne[i][j];
        }
      }
    } 
  } 

  if (debug == 1) {
    printMas( temp, xPoints, yPoints );
  }

  if (gnu == 1) { 
    long iterations = calcTime / dt;
    long passedIteration;
    double showTime = 5.0;

    if (iterations > (showTime / dt) ) {
      passedIteration = iterations / (showTime / dt);
    } else {
      passedIteration = 1;
    }
     
    if ( debug == 1) {
      printf ("\n iterations = %ld, passed = %ld, %ld ",
        iterations, passedIteration, (long)(showTime / dt));
    } 

    if ( calcIterat % passedIteration == 0)
    {
      makeFinalFileSETUP_phase0(temp);
      makeFinalFileDATA_phase1(temp);
    } 
  }
}

void makeArrayForSendToGNU(double* masFrom, double** masTo, int total) {
  
  printf("\nHELLO!!!!!\n");
  long size = yPoints + total*2; 
  long cX = 0, cY = 0; 
  long shift = 0;
  long i, j;
  for (i = 0; i < size; i++) {
    for (j = 0; j < xPoints; j++) {
      
      if (debug == 1 && flagForPrintArrayTOGNUANDSOLV) {
        printf("\ni = %ld,%ld cX = %ld, shify = %ld, size= %ld", 
          i, j, cX, shift, size);
      }

      if ((i == 0) || (i == size-1)) {
        //
      } else 
      if ( shift == 0) {
        masTo[cX][j] = masFrom[i*xPoints +j]; 
        printf("++++++++++");
      } else {
        if (shift == 2) {
          //
        } else {
          //
        } 
      }
    }

    if ( shift == 0) {
      if (i != 0) {
        if ((cX < yPoints) && (shift == 0)) cX ++;
        if ( cX % (yPoints/total) == 0 ) {
          shift = 2;
        }
      }
       
    } else { 
      shift --;
    }
  }   

  if (debug == 1) {
    printf("\n Make Array For GNU! \n");
    for (i = 0; i < yPoints; i++) {
      for (j = 0; j < xPoints; j++) { 
        printf("%9.4f ", masTo[i][j]);
      }
      printf("\n");
    }
  }
}

void makeArrayForSendToSolve(double* masTo, double** masFrom, int total) {
  
  long cX = 0, cY = 0; 
  long shift = 0;
  long i, j;
  for (i = 0; i < yPointsMPI; i++) {
    for (j = 0; j < xPoints; j++) {
       
      if (debug == 1 && flagForPrintArrayTOGNUANDSOLV) {
        printf("\ni = %ld,%ld cX = %ld, shify = %ld, size= %ld", 
          i, j, cX, shift, yPointsMPI);
      }

      if ((i == 0) || (i == yPointsMPI-1)) {
        masTo[i*xPoints +j] = 0.0;
      } else 

      if ( shift == 0) {
        masTo[i*xPoints +j] = masFrom[cX][j]; 
      } else {
        if (shift == 2) {
          masTo[i*xPoints +j] = masFrom[cX][j];
        } else {
          masTo[i*xPoints +j] = masFrom[cX-1][j];
        } 
      }
    }

    if ( shift == 0) {
      if (i != 0) {
        if ((cX < yPoints) && (shift == 0)) cX ++;
        if ( cX % (yPoints/total) == 0 ) {
          shift = 2;
        }
      }
       
    } else { 
      shift --;
    }
  } 

  if (debug == 1) {
    printf("\n Make Array For solve! \n");
    for (i = 0; i < yPointsMPI; i++) {
      for (j = 0; j < xPoints; j++) { 
        printf("%6.2f ", masTo[i*xPoints +j]);
      }
      printf("\n");
    }
  }
}

void findMinMaxGnuPlot(double** temp){

  //printf("\n FIND MAX MIN TEMP"); 
  long i, j; 

  maxT = -100000.0;
  minT = 100000.0;

  for (i = 0; i < yPoints; i++) {
    for (j = 0; j < xPoints; j++) { 
      if (temp[i][j] > maxT) {
        maxT = temp[i][j];
      }
      if (temp[i][j] < minT) {
        minT = temp[i][j];
      } 
    } 
  }
   
  //printf("\n max = %f, min = %f \n", maxT, minT);
}

void makeFinalFileFileOpen() {
  // calculate pDx, pDy
  if (numPoints > 1000) {
    printf("\n NumPoints more than 1000 \n");

    int x = sqrt( 1000 / ((length) * (hight)) );
    printf("\n x= %d", x);

    printDx = (double)(length / ( ((length) * x) -1.0) );
    printDy = (double)(hight / ( ((hight) * x) -1.0) );
    printf("\n dx= %f(%f), dy= %f(%f) ", 
      (double)printDx, (double)printDx*(length*((int)x) -1.0), 
      (double)printDy, (double)printDy*(hight*((int)x) -1.0) );

    printXPoints = length / printDx;
    printYPoints = hight / printDy;

    printf("\n printfDx= %ld, printDy= %ld", printXPoints, printYPoints);

  } else {
    printf("\n NumPoints less than 1000 \n");
    printXPoints = xPoints;
    printYPoints = yPoints;

    printDx = dx;
    printDy = dy;
  }
 
  printf("\n OPEN FILE \n");
  file = fopen("file.txt", "w");
  fprintf(file, "reset \n");
  fprintf(file, "clear \n");
  fprintf(file, "plot 2 with lines \n");
  fprintf(file, "set size ratio %f \n", (double)(hight / length) );
  fprintf(file, "set nokey \n");
  fprintf(file, "set xrange [-1:%ld] \n", (long)(length/printDx));
  fprintf(file, "set yrange [-1:%ld] \n", (long)(hight/printDy));
}

void makeFinalFileSETUP_phase0(double** temp) {
  findMinMaxGnuPlot(temp);
  fprintf(file, "set cbrange [%f:%f] \n", minT -1, maxT +1);
  fprintf(file, "set palette defined (%ld \"#0000FF\", %ld \"#FF0000\") \n", 
    (long)minT -1, (long)maxT +1 );
  fprintf(file, "plot \"-\" using 2:1:3 with image \n");
}
 
void makeFinalFileDATA_phase1(double** temp) {
  long i, j; 

  for (i = 0; i < printYPoints; i++) {
    for (j = 0; j < printXPoints; j++) { 
      
      fprintf(file, "%ld %ld %f \n", i, j, temp[(int)(i*(printDy/dy))][(int)(j*(printDx/dx))]);

    } 
  }
  fprintf(file, "EFO \n pause 0.1 \n"); 
}

void makeFinalFileEND_phase2() {
  printf("\n FILE CLOSED");
  fclose(file);
}

void inicialBorderConditionsOne(double** border, double temp) {
  long i, j;

  for (i = 0; i < yPoints; i++) {
    border[i][0] = temp;
  }

  //border[(long)yPoints/2][(long)xPoints/2] = temp;
}

void inicialBorderConditionsTwo(double** border, double temp) {
  long i, j;

  for (i = 0; i < yPoints; i++) {
    border[i][(long)xPoints-1] = temp;
  } 
}

void applyBorderOneToTemp(double** borderOne, double** temp) {
  long i, j; 

  for (i = 0; i < yPoints; i++) {
    for (j = 0; j < xPoints; j++) {
      if (borderOne[i][j] != defCondOne) {
        temp[i][j] = borderOne[i][j];
      }
    } 
  }
}

void inicialBeginTemp(double** x, double temp) {
  long i, j; 

  for (i = 0; i < yPoints; i++) {
    for (j = 0; j < xPoints; j++) {
      x[i][j] = temp;
    } 
  }
}

double** make2DArray(long arraySizeX, long arraySizeY) {
  double** theArray;
  theArray = (double**) malloc(arraySizeY*sizeof(double*));
  long i;
  for (i = 0; i < arraySizeY; i++) {
    theArray[i] = (double*) malloc(arraySizeX*sizeof(double));
  }
  return theArray;
}

void compareTwoArray(double** x, double** y, long ix, long jy){
  long flag = 0;
  long i, j;
  for (i = 0; i < jy; i++) {
    for (j = 0; j < ix; j++) {
      if (x[i][j] != y[i][j]) {
        flag ++;
      }
    }
  }
  if (flag == 0) {
    if (debug == 1) {
      printf("\nCompare was similar ! +++++++++++\n");
    }
  } else {
    printf("\nCompare wasn't similar ! ------------- [%ld]\n", flag);
    if (debug == 1) {
      printMas(x, ix, jy);
      printMas(y, ix, jy);
      printf("\nPRINTED SOLUTION X and Y\n");
    }
  }
}

void printMas(double** x, long ix, long jy){
  long i, j;

  printf("\n print massive i=%ld, j=%ld. %ld \n", ix, jy, sizeof(x));
  for (i = 0; i < jy; i++) {
    for (j = 0; j < ix; j++) {
      printf("%9.4f ", x[i][j]);
    }
    printf("\n");
  }
  printf(" end print massive \n");
}

void printMasI(double* x, long ix, long jy){
  long i, j;

  printf("\n print massive i=%ld, j=%ld. \n", ix, jy);
  for (i = 0; i < jy; i++) {
    for (j = 0; j < ix; j++) {
      printf("%9.4f ", x[i*ix +j]);
    }
    printf("\n");
  }
  printf(" end print massive \n");
}

 
