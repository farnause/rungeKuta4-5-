int rk45(double *ev, double *e, double *t, double *x,  double *h, double hmin, double hmax, double tol, int n, int (*camp)(int n, double t, double *x, double *f, double *prm), double *prm);

int flux (double *ev, double *e, double *t, double *x, double *h, double T, double hmin,double hmax, double tol, int npasmax, int n,int (*camp)(int n, double t, double *x, double *f, double *prm),double *prm);