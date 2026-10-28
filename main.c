
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Funcio.h"
#include "rkf45.h"
#include "otros.h"



int main(){


    double t,hmin,hmax,tol,h,A,x0,y0;
    int n=2;
    double x[n];
    double e=0;
    double ev[n];

    Start(ev,n);

    x0=1;
    y0=2;

    x[0]=x0;
    x[1]=y0;

    double a=0;
    double T=1;

    t=0;//valor inicial de t
    tol=1e-10; //tolerancia que se le permite al error

    //valores entre los que h se puede encontrar
    hmin=tol;
    hmax=0.1;
    
    h=0.01;//paso inicial

    A=rk45(ev,&e,&t, x,&h,hmin, hmax,tol,n,camp,&a);

    if(A==-1){
        printf("Ha ocurrido un error\n");
    }else if (A==-2){
        printf("El valor de h se ha reducido por debajo de hmin. T final: %.16G\n", t);
    }

    printf("Posicion final calculada: (%.16G, %.16G), t: %.16G\n", x[0], x[1],t);
    printf("Posicion final real: (%.16G, %.16G)", exp(a*t)*(x0*cos(t)-y0*sin(t)), exp(a*t)*(x0*sin(t)+y0*cos(t)));


    return 0;
}
