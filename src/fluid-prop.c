#include <stdio.h>
#include "libFluidPropC.h"

int main(){
	double T,P,rho,alpha,comp,lambda,kappa,nu,cp,psi,lewis;

	T = 300.;
	P = 1.0e5;
	getCoolProp("Water",T,P,&rho,&alpha,&comp,&lambda,&kappa,&nu,&cp,&psi,&lewis);

	printf("T: %lf K\tP: %lf bar\n",T,P*1e-5);
	printf("mass density       (kg/m^3): %.3g \n",rho);
	printf("expansion coefficient (1/K): %.3g \n",alpha);
	printf("heat conductivity   (W/m K): %.3g \n",lambda);

	return 0;
}
