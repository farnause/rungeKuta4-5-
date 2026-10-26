
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Funcio.h"
#include "rkf45.h"
#include "otros.h"



int main(){

    //HE AFEGIT LA OPCI� DE QUE MOSTRI COTA ERROR MODUL I VECTOR ERROR

    double t,hmin,hmax,tol,h,A;
    int n=2;
    double x[n];
    double e=0;
    double ev[n];

    Start(ev,n);

    x[0]=1;//x0
    x[1]=2;//y0

    double a=1;
    double T=1;

    t=0;
    tol=1e-10;
    hmin=tol;
    hmax=0.1;
    h=0.01;

    A=rk45(ev,&e,&t, x,&h,hmin, hmax,tol,n,camp,&a);

    if(A==-1){
        printf("Ha ocurregut un error\n");
    }else if (A==-2){
        printf("S'ha escurçat el pas per sota de hmin. T final: %.16G\n", t);
    }

    printf("Posicio final: (%.16G,%.16G)", x[0], x[1]);


    return 0;
}
