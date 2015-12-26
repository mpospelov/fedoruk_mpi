#include <stdio.h>  
#include <string.h> 
#include <stdlib.h> 

#define debug 1 // 1 - on, 0 - off
#define gnu 1 // 1 - on, 0 - off

#define defCondOne -999999 // inicialize first border condition

double calcTime;
double dx, dy, dt;
double length, hight;
double xPoints, yPoints;
double at;

double** tempMas;
double** tempMasPrev;

double** borderConditionOne;
 
double** make2DArray(long arraySizeX, long arraySizeY);

void inicialBeginTemp(double **, double);

void printMas(double**, long, long);
long R(long , long , long );
void solve(double** , double** , double** );
void inicialBorderConditions( double** , double );
void applyBorderOneToTemp(double** ,double**);

// gcc -g -o 4 4.c

int main(){

	calcTime = 8800; // seconds

	dx = 0.01;
	dy = 0.01;
	dt = 0.1;

	length = 2.0;
	hight = 1.0;

	at = 420.0 / (235.0 * 10500.0);

	xPoints = 1 + length/dx;
	yPoints = 1 + hight/dy;

	tempMas = make2DArray(xPoints, yPoints);
	tempMasPrev = make2DArray(xPoints, yPoints); 

	inicialBeginTemp(tempMas, 100.0);
	inicialBeginTemp(tempMasPrev, 100.0);

	borderConditionOne = make2DArray(xPoints, yPoints); 

	inicialBeginTemp(borderConditionOne, defCondOne);

	inicialBorderConditions(borderConditionOne, 0.0); // temp in x[i][0] = temp

	printMas( tempMas, xPoints, yPoints );
	printMas( tempMasPrev, xPoints, yPoints );
	printMas( borderConditionOne, xPoints, yPoints );

	printf("\n\n START SOLVE \n\n");

	applyBorderOneToTemp(borderConditionOne, tempMasPrev);

	long i;
	int flag = 0;
	for (i = 0; i <= (calcTime / dt); i++)
	{
		if (flag == 0) {
			solve(tempMas, tempMasPrev, borderConditionOne); 
			flag = 1;
		} else
		{
			solve(tempMasPrev, tempMas, borderConditionOne);
			flag = 0;
		} 
	}

 

	return 0;
}

long R(long i, long k, long max) {
	long ans = i + k;
	//printf("max=%ld\n", max);
	if ((ans < max) && (ans > 0)) {
		//printf("ans=%ld\n",ans);
		return ans;
	
	} else {

		if (ans >= max) {
			//printf("ans=%ld\n",ans);
			return max - 1;
		} else {
			//printf("ans=%ld\n",ans);
			return 0;
		}
	}
}

void solve(double** temp, double** tempPrev, double** borderOne) {

	long i, j, k;
	printf("$\n");
	for (i = 0; i < yPoints; i++) {
		for (j = 0; j < xPoints; j++){

			if (borderOne[i][j] == defCondOne) { // without first border condition

				double tdx, tdy;
				
				tdy = (double) ( tempPrev[i][ R(j, 1, xPoints)] 
					- 2.0 * tempPrev[i][j] 
					+ tempPrev[i][ R(j, -1, xPoints)] ) / (double) ( dx * dx );

				tdx = (double) ( tempPrev[ R(i, 1, yPoints)][j] 
					- 2.0 * tempPrev[i][j] 
					+ tempPrev[ R(i, -1, yPoints)][j] ) / (double) ( dy * dy );


				temp[i][j] = (double) dt * at * ( tdx + tdy ) + tempPrev[i][j];

			} else { // first border condition
				
				temp[i][j] = borderOne[i][j];
			}
		} 
	} 

	printMas( temp, xPoints, yPoints );
}

void inicialBorderConditions(double** border, double temp) {
	long i, j;

	for (i = 0; i < yPoints; i++) {
		border[i][0] = temp;
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

 