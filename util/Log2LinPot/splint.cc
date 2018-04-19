#include <vector>
#include <iostream>
#include <stdlib.h>
using namespace std;

void splint(double* xa,double* ya, double* y2a, const int n, double x, double *y)
{
    int klo=0;
    int khi=n-1;
    
    while (khi-klo > 1){
        int k=(khi+klo)>>1;
        if(xa[k]>x)khi=k;
        else klo=k;
    }
    double h=xa[khi]-xa[klo];
    if(h==0.){
        cout<<"Bad xa input"<<endl;
        exit(0);
    }
    double a=(xa[khi]-x)/h;
    double b=(x-xa[klo])/h;
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.;
}

