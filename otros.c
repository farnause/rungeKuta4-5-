
#include <stdio.h>
#include <stdlib.h>

double Abs(double x){
    if(x>=0){
        return x;
    }else{
        return -x;
    }
}
double minomax(double x, double y, int a){//si a=-1 et torna el m�nim, si a=1 et torna el m�xim

    if(x<y && a==-1){
         return x;
    }else if (x>y && a==-1){
        return y;
    } if(x<y && a==1){
         return y;
    }else if (x>y && a==1)
        return x;
    return x;//cas que x=y
}
void Start (double *v, int n){
    for (int i=0; i<n; i++){
        v[i]=0;
    }
}
