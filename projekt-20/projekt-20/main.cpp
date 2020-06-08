#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rk4.h"

#define g 9.81           //m/s^2
#define ro_powietrza 1.2 //kg/m^3 
#define Pi 4*atan(1.)
#define c 0.47   //wspolczynnik oporu powietrza dla kuli

double S;    //powierzchnia na ktora dziala sila oporu
double m;

void rhs_fun1(double t, double* X, double* A);
void rhs_fun2(double t, double* Y, double* B); 

void main()
{
	double V0;  //m/s
	double V_powietrza; //m/s, rosnie liniowo wraz z wysokoscia
	double alfa;    
	double Vx0;
	double Vy0;
	double Vx, Vy; 
	double x0 = 0.0;
	double y0 = 0.0;
	int n = 2;
	double h = 0.1;

	printf("Podaj wartosc masy kuli (w kg), wartosc powierzchni na jaka dziala sila oporu (w m^2): \n");
	scanf("%lf %lf", &m, &S);
	printf("Podaj wartosci predkosci poczatkowej kuli (w m/s)  i wartosc kata pod jakim kula zostala wystrzelona (w stopniach): \n");
	scanf("%lf %lf", &V0, &alfa);


	alfa = alfa*Pi/180.;    
	Vx0 = V0 * cos(alfa);
	Vy0 = V0 * sin(alfa);
	Vx = Vx0 - V_powietrza;  //predkosc kuli wzgledem wiatru
	V_powietrza = y0*5.; 

	double* X = (double*)malloc(n*sizeof(double)); 
	double* X1 = (double*)malloc(n*sizeof(double));   
	double* Y = (double*)malloc(n*sizeof(double));    
	double* Y1 = (double*)malloc(n*sizeof(double)); 
    X[0] = Vx;   
	X[1] = x0;
	Y[0] = Vy;
	Y[1] = y0;


	FILE* f;
	f = fopen("dane.txt", "w");

	for(double t=0;y >= 0;t+=h)
	{
		vrk4(t, X, h, n, rhs_fun1, X1);
		X[0] = X[1];
		X[1] = X1[1];
	
		vrk4(t, Y, h, n, rhs_fun2, Y1);
		Y[0] = Y[1];
		Y[1] = Y1[1];

		fprintf(f, "%lf\t%lf\t%lf\t%lf\t%lf\n", t, X1[0], X1[1], Y1[0], Y1[1]);
	}

	fclose(f);
	free(X);
	free(X1);
	free(Y);
	free(Y1);

}

void rhs_fun1(double t, double* X, double* A) 
{
	A[0] = -(ro_powietrza*X[0]*X[0]*S*c)/(2*m);
	A[1] = X[0];						
}


void rhs_fun2(double t, double* Y, double* B) 
{
	B[0] = ((ro_powietrza*Y[0]*Y[0]*S*c)/(2*m))-g;
	B[1] = Y[0];						
}

