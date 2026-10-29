
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Funcio.h"
#include "rkf45.h"
#include "otros.h"


int main(){


    int n=2;
    int npasmax;
    double x[n];
    double e=0;
    double ev[n];
    double a,t,T,tol,hmax,hmin,h,x0,y0;


    Start(ev,n);

    //condiciones iniciales
    x0=1;
    y0=2;
    t=0;//tiempo inicial

    a=1;//valor del parámetro
    T=3;//tiempo final
    tol=1E-10;//tolerancia que se le permite al error
    h=0.01;//valor inicial de h

    //valores entre los que h se puede encontrar
    hmin=tol;
    hmax=0.1;

    npasmax=1E6; //número de pasos máximos

    x[0]=x0;
    x[1]=y0;
    
    if(flux(ev,&e, &t,x,&h,T,hmin,hmax,tol,npasmax,n,camp,&a)<0){
        printf("No se ha podido integrar numericament hasta el final\n");
        printf("Posicion final calculada: (%.16G, %.16G), t: %.16G\n", x[0], x[1]);
        printf("Posicion final real: (%.16G, %.16G)", exp(a*t)*(x0*cos(t)-y0*sin(t)), exp(a*t)*(x0*sin(t)+y0*cos(t)));
    }else{
        printf("Posicion final calculada: (%.16G, %.16G), t: %.16G\n", x[0], x[1]);
        printf("Posicion final real: (%.16G, %.16G)", exp(a*t)*(x0*cos(t)-y0*sin(t)), exp(a*t)*(x0*sin(t)+y0*cos(t)));
    }

    


    return 0;
}
