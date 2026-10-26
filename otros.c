
#include <stdio.h>
#include <stdlib.h>

double Abs(double x){//función propia de valor absoluto
    if(x>=0){
        return x;
    }else{
        return -x;
    }
}
double minomax(double x, double y, int a){//si a=-1 la función devuelve el mínimo, si a=1 te devuelove el máximo

    if(x<y && a==-1){
         return x;
    }else if (x>y && a==-1){
        return y;
    } if(x<y && a==1){
         return y;
    }else if (x>y && a==1)
        return x;
    return x;
}
void Start (double *v, int n){//sirve para "inicializar" un vector con sus elementos en 0
    for (int i=0; i<n; i++){
        v[i]=0;
    }
}
