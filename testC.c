#include <stdio.h>
#include "libFluidPropC.h"

int main(){
	double T,P,rho,alpha,comp,lambda,kappa,nu,cp;

	T = 300.;
	P = 1.0e5;
	getCoolProp("Water",T,P,&rho,&alpha,&comp,&lambda,&kappa,&nu,&cp);
	printf("T: %lf\tP: %lf\n",T,P);
	printf("rho: %.3lf\n",rho);

	return 0;
}
