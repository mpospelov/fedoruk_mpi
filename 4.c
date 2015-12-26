#include <stdio.h>  
#include <string.h> 
#include <stdlib.h> 
#include <math.h>

#define debug 0 // 1 - on, 0 - off
#define gnu 1 // 1 - on, 0 - off

#define defCondOne -999999 // inicialize first border condition

double calcTime;
double dx, dy, dt;
double length, hight;
double xPoints, yPoints;
double at;

double printDx, printDy;
long printXPoints, printYPoints;

long numPoints;

double maxT = -100000.0;
double minT = 100000.0;

FILE *file;

double** tempMas;
double** tempMasPrev;

double** borderConditionOne;
double** borderConditionTwo;
 
double** make2DArray(long arraySizeX, long arraySizeY);

void inicialBeginTemp(double **, double);

void printMas(double**, long, long);
long R(long, long, long, long );
void solve(double** , double** , double** , double** , long);
void inicialBorderConditionsOne( double** , double );
void inicialBorderConditionsTwo( double** , double );
void applyBorderOneToTemp(double** ,double**);

void findMinMaxGnuPlot(double** );
void makeFinalFileFileOpen();
void makeFinalFileSETUP_phase0(double** );
void makeFinalFileDATA_phase1(double** );
void makeFinalFileEND_phase2();

// gcc -g -o 4 4.c -lm
// ./4 10000 512

int main( int argc, char **argv ){

	dt = 0.1;

	length = 2.0;
	hight = 1.0;

	if (argc < 3)
	{
		printf("\nError input parameters!\n");
		 	
		calcTime = 300; // seconds

		dx = 0.1;
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

	inicialBeginTemp(tempMas, 100.0);
	inicialBeginTemp(tempMasPrev, 100.0);

	borderConditionOne = make2DArray(xPoints, yPoints); 
	borderConditionTwo = make2DArray(xPoints, yPoints); 

	inicialBeginTemp(borderConditionOne, defCondOne);
	inicialBeginTemp(borderConditionTwo, defCondOne);

	inicialBorderConditionsOne(borderConditionOne, 0.0); // temp in x[i][0] = temp
	inicialBorderConditionsTwo(borderConditionTwo, 100.0); // heat temp in x[i][last] = temp

	if (debug == 1) {
		printMas( tempMas, xPoints, yPoints );
		printMas( tempMasPrev, xPoints, yPoints );
		printMas( borderConditionOne, xPoints, yPoints );	
	} 

	printf("\n\n START SOLVE \n\n");

	applyBorderOneToTemp(borderConditionOne, tempMasPrev);

	if (gnu == 1) {
		findMinMaxGnuPlot(tempMasPrev);	
		makeFinalFileFileOpen();
	}

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

					double tdx, tdy;
					
					tdy = (double) ( tempPrev[i][ R(j, 1, 0, xPoints)] 
						- 2.0 * tempPrev[i][j] 
						+ tempPrev[i][ R(j, -1, 0, xPoints)] ) / (double) ( dx * dx );

					tdx = (double) ( tempPrev[ R(i, 1, 0, yPoints)][j] 
						- 2.0 * tempPrev[i][j] 
						+ tempPrev[ R(i, -1, 0, yPoints)][j] ) / (double) ( dy * dy );


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

	border[(long)yPoints/2][(long)xPoints/2] = temp;
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

void printMas(double** x, long ix, long jy){
	long i, j;

	printf("\n print massive i=%ld, j=%ld. %ld \n", ix, jy, sizeof(x));
	for (i = 0; i < jy; i++) {
		for (j = 0; j < ix; j++) {
			printf("%f ", x[i][j]);
		}
		printf("\n");
	}
	printf(" end print massive \n");
}

 