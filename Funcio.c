#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int camp (int n, double t, double *x, double *f, double *prm){
    if(n!=2){
        printf("La funcio esta pensada per n=2!!!");
        return -1;
    }else if(prm==NULL){
        f[0]=x[0]-x[1];
        f[1]=x[0]+x[1];
    }else{
        f[0]=(*prm)*x[0]-x[1];
        f[1]=x[0]+(*prm)*x[1];
    } 

    return 1;
}
