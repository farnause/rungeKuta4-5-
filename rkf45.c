#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Funcio.h"
#include "otros.h"


int rk45(double *ev, double *e, double *t, double *x,  double *h, double hmin, double hmax, double tol, int n, int (*camp)(int n, double t, double *x, double *f, double *prm), double *prm)
{
	/*Coeficients RK45 (estatics)*/
	static double c[6] = {0., 1./4., 3./8., 12./13., 1., 1./2.};
	static double A[15] = {1./4., 3./32., 9./32., 1932./2197.,
		-7200./2197.,7296./2197., 439./216., -8., 3680./513.,
		-845./4104., -8./27., 2., -3544./2565., 1859./4104., -11./40.};
	static double b[6] = {25./216., 0., 1408./2565., 2197./4104., -1./5.,0.};
	static double bp[6] = {16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.};
	static double E[6] = {-1./360., 0., 128./4275., 2197./75240., -1./50., -2./55.};
	double p=0.9;
	double eps=0;//epsilon per calcular distancia entre x7 i x8
	double aux; //integer que si no no m'accepta el �ndex de A en el calcul de y2

	*(h)=(Abs(*(h))/(*(h)))*minomax(hmax,minomax(hmin,Abs(*(h)),1),-1);//Per seguretat

	double *f=malloc(n* sizeof(double));//entenc que dim arribada=dim sortida-1 (contem la variable t) de f,(Dim camp =dim x +1 de t)??
	double *y1=malloc(n* sizeof(double));//vecotor auxiliar per calcular lo que se li afegeix a x en el calcul de cada kapa
	double *y2=malloc(n* sizeof(double));//vecotor auxiliar per calcular beta*kapa anterior
	double *y3=malloc(n* sizeof(double));//per calcul error
	double *y=malloc(n* sizeof(double));//vector auxiliar en el que s'emmmagatzema x+y1
	double *K=malloc((n*6)*sizeof(double));//6 etapes
	double *x4=malloc(n* sizeof(double));//entenc que dim arribada=dim sortida-1 (contem la variable t) de f,(Dim camp =dim x +1 de t)??
	double *x5=malloc(n* sizeof(double));
	double *epsv=malloc(n* sizeof(double));

    camp(n,*t,x,f,prm);
    for(int i=0;i<n;i++){
        K[i]=f[i]; //K ho he pensat com una matriu(o un vector) en que cada fila �s una kapa
    }//"K[0]=f(t,x)"
	while(1){
        for(int l=1; l<6; l++){//Aqu� s'ha de canviar el n�mero de etapes en cas de voler canviar de rungekuta
            aux=(l*l-l)/2.;
            Start(y1,n);
            for(int s=0; s<l;s++){//sumatori
                for(int g=0; g<n; g++){//calculo beta*kapa anteriro
                    y2[g]=(*h)*A[(int)aux+s]*K[n*s+g];
                    y1[g]=y1[g]+y2[g];
                }
            }
            for (int j=0; j<n; j++){
                y[j]=x[j]+y1[j];
            }
            camp(n,*t+c[l],y,f,prm);

            for(int i=0;i<n;i++){
                K[l*n+i]=f[i]; //K ho he pensat com una matriu(o un vector) en que cada fila �s una kapa
            }
            //printf("\n");
        }
        //C�LCUL x4 i x5
        //Puc reutilitzar y1, y2 i y que tenen mida 6
        Start(x4,n);
        Start(x5,n);
        Start(epsv,n);


        for(int s=0; s<6;s++){
            for(int i=0; i<n; i++){
                y2[i]=(*h)*b[s]*K[n*s+i];
                y1[i]=(*h)*bp[s]*K[n*s+i];
                y3[i]=(*h)*E[s]*K[n*s+i];
                x4[i]=x4[i]+y2[i];
                x5[i]=x5[i]+y1[i];
                epsv[i]=epsv[i]+y3[i];//calculem el error.
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
            if(eps<(tol/(2*2*2*2*2))){//acotem el n�mero pel cual multiplequem |h| per 2, al calcular el nou h.
                *(h)=(Abs(*(h))/(*(h)))*minomax(hmax,2*p*Abs(*(h)),-1);
                return 1;
            }
            *(h)=(Abs(*(h))/(*(h)))*minomax(hmax,minomax(hmin,Abs(*(h))*p*pow(tol/eps,1./5.),1),-1);//POW POTSER ES MASSA LENT? AIX� POW ESTA B�?
            return 1;

        }
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
            return;//En teoria s'gauria de tenir en compte per canviar hmin
        }

        printf("%G*%G*(%G/%G)^1/5, ",Abs(*(h)),p,tol,eps);
        *(h)=(Abs(*(h))/(*(h)))*minomax(hmin,Abs(*(h))*p*pow(tol/eps,1./5.),1);
        printf("%G=h\n",*(h));
        eps=0;
	}
	return 0;
}


