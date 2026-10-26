#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Funcio.h"
#include "otros.h"


int rk45(double *ev, double *e, double *t, double *x,  double *h, double hmin, double hmax, double tol, int n, int (*camp)(int n, double t, double *x, double *f, double *prm), double *prm)
{
	/*Coeficientes RK45 (estaticos)*/
	static double c[6] = {0., 1./4., 3./8., 12./13., 1., 1./2.};
	static double A[15] = {1./4., 3./32., 9./32., 1932./2197.,
		-7200./2197.,7296./2197., 439./216., -8., 3680./513.,
		-845./4104., -8./27., 2., -3544./2565., 1859./4104., -11./40.};
	static double b[6] = {25./216., 0., 1408./2565., 2197./4104., -1./5.,0.};
	static double bp[6] = {16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.};
	static double E[6] = {-1./360., 0., 128./4275., 2197./75240., -1./50., -2./55.};
	double p=0.9;
	double eps=0;//varaible para calcular la distancia entre x7 i x8
	double aux; //varaible auxiliar en el calculo de y2

	*(h)=(Abs(*(h))/(*(h)))*minomax(hmax,minomax(hmin,Abs(*(h)),1),-1);//Para asegurarse que h esta eb (hmin, hmax)

	double *f=malloc(n* sizeof(double));
	double *y1=malloc(n* sizeof(double));//vecotor auxiliar para calcular lo que se le añade a x
	double *y2=malloc(n* sizeof(double));//vecotor auxiliar para calcular beta*kapa anterior
	double *y3=malloc(n* sizeof(double));//por el calculo del error
	double *y=malloc(n* sizeof(double));//vector auxiliar en el que se alamcena x+y1
	double *K=malloc((n*6)*sizeof(double));//aquí es necesario saber el número de etapas del rkf (en este caso 6)
	double *x4=malloc(n* sizeof(double));//aproximación
	double *x5=malloc(n* sizeof(double));//valor de mayor "precisión" con el que se calcula el error
	double *epsv=malloc(n* sizeof(double));

    camp(n,*t,x,f,prm);
    for(int i=0;i<n;i++){
        K[i]=f[i]; //K respresenta una matriz en que cada hilera es una etapa
    }//"K[0]=f(t,x)"
	while(1){

        //CALCULO DE LAS ETAPAS
        for(int l=1; l<6; l++){//Aquí tambien se ha de tener en cuenta el número de etapas (6)
            aux=(l*l-l)/2.;
            Start(y1,n);
            for(int s=0; s<l;s++){
                for(int g=0; g<n; g++){
                    y2[g]=(*h)*A[(int)aux+s]*K[n*s+g];
                    y1[g]=y1[g]+y2[g];
                }
            }
            for (int j=0; j<n; j++){
                y[j]=x[j]+y1[j];
            }
            camp(n,*t+c[l],y,f,prm);

            for(int i=0;i<n;i++){
                K[l*n+i]=f[i]; 
            }
        }

        //CÁLCULO x4 y x5
        Start(x4,n);
        Start(x5,n);
        Start(epsv,n);

        for(int s=0; s<6;s++){//Aquí tambien se ha de tener en cuenta el número de etapas (6)
            for(int i=0; i<n; i++){
                y2[i]=(*h)*b[s]*K[n*s+i];
                y1[i]=(*h)*bp[s]*K[n*s+i];
                y3[i]=(*h)*E[s]*K[n*s+i];
                x4[i]=x4[i]+y2[i];
                x5[i]=x5[i]+y1[i];
                epsv[i]=epsv[i]+y3[i];//cálculo error
                if(isnan(x5[i])||isnan(x5[i])){
                    free(f);
                    free(y1);
                    free(y2);
                    free(y3);
                    free(y);
                    free(K);
                    free(x4);
                    free(x5);
                    free(epsv);
                    return -1;
                }
            }
        }

        //cálculo modulo error
        for(int i=0; i<n; i++){
            eps=eps+Abs(epsv[i]);
            ev[i]=Abs(epsv[i])+ev[i];
        }
        *e=*e+eps;

        if(eps<tol){
            *(t)=*(h)+*(t);
            for(int i=0; i<n; i++){
                x[i]=x5[i]+x[i];
            }
            free(f);
            free(y1);
            free(y2);
            free(y3);
            free(y);
            free(K);
            free(x4);
            free(x5);
            free(epsv);
            if(eps<(tol/(2*2*2*2*2))){//acotamos el número por el cual se multiplica |h| por 2, al calcular el nuevo h.
                *(h)=(Abs(*(h))/(*(h)))*minomax(hmax,2*p*Abs(*(h)),-1);
                return 1;
            }
            *(h)=(Abs(*(h))/(*(h)))*minomax(hmax,minomax(hmin,Abs(*(h))*p*pow(tol/eps,1./5.),1),-1);
            return 1;

        }

        //comprovación del h obtenido
        if(Abs(*(h))<=hmin){
            *(t)=*(h)+*(t);
            for(int i=0; i<n; i++){
                x[i]=x5[i]+x[i];
            }
            *(h)=hmin;
            free(f);
            free(y1);
            free(y2);
            free(y3);
            free(y);
            free(K);
            free(x4);
            free(x5);
            free(epsv); 
            return -2; //Mostramos error si h llega a ser menor que hmin
        }

        *(h)=(Abs(*(h))/(*(h)))*minomax(hmin,Abs(*(h))*p*pow(tol/eps,1./5.),1);//paso de seguridad para el nuevo h
        eps=0;
	}
	return 0;
}


