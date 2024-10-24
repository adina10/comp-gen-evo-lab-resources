#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double sine(double x);
double powerd(double x, int y);
double factorial_div(double value, int x);

int main(void){

int TIME, PERIOD;
long double WIDTH;

FILE * out, * par;

par = fopen("parst", "r"); //input file

fscanf(par, "%d", &TIME); // number of generations
fscanf(par, "%d", &PERIOD); // number of generations to go around the "circle"
fscanf(par, "%Lf", &WIDTH); // smax

char filename[64];
sprintf(filename, "trajectory.%dx%.2Lf", PERIOD, WIDTH);

out = fopen(filename, "w"); // output file

WIDTH *=2.0;

long double freqs;
double y, freq;
int k;

for ( k = 0; k < TIME; k++){
y = k  * 2 * 3.1415926 / PERIOD;     
    freq = 0.0 + (WIDTH/2.00) * sine(y); // the sine function
    if(k < PERIOD)  printf("%f \n", freq); // just to see one oscilation on the screen
    fprintf(out, "%f \n" , freq); // saves trajectory
}

return 0;
}

/*----------------------------------------------------------------------------*/
double sine( double x)
{
int i=0;
int j=1;
int sign=1;
double y1 = 0.0;
double diff = 1000.0;

if (x < 0.0){
x = -1 * x;
sign = -1;
}

while ( x > 360.0*3.1415926/180){
x = x - 360*3.1415926/180;
}

if( x > (270.0 * 3.1415926 / 180) ){
sign = sign * -1;
x = 360.0*3.1415926/180 - x;
}

else if ( x > (180.0 * 3.1415926 / 180) ){
sign = sign * -1;
x = x - 180.0 *3.1415926 / 180;
}

else if ( x > (90.0 * 3.1415926 / 180) ){
x = 180.0 *3.1415926 / 180 - x;
}

while( powerd( diff, 2) > 1.0E-16 ){
i++;
diff = j * factorial_div( powerd( x, (2*i -1)) ,(2*i -1));
y1 = y1 + diff;
j = -1 * j;
}

return ( sign * y1 );
}

/*----------------------------------------------------------------------------*/
double powerd( double x, int y){
int i=0;
double ans=1.0;

if(y==0) return 1.000;
else{
while( i < y){
i++;
ans = ans * x;
}
}
return ans;
}

/*----------------------------------------------------------------------------*/
double factorial_div( double value, int x){
if(x == 0) return 1;
else{
while( x > 1){
value = value / x--;
}
}
return value;

}