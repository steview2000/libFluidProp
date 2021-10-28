#include <stdio.h>
#include "libFluidPropC.h"

int main(){
	double T,P,rho,alpha,comp,lambda,kappa,nu,cp,psi,lewis;
	char flag[100];

	sprintf(flag,"INCOMP");
	T = 300.;
	P = 1.0e5;
	getCoolProp("MGL[0.50]",T,P,&rho,&alpha,&comp,&lambda,&kappa,&nu,&cp,&psi,&lewis,flag);

	printf("T: %lf K\tP: %lf bar\n",T,P*1e-5);
	printf("mass density       (kg/m^3): %.3g \n",rho);
	printf("expansion coefficient (1/K): %.3g \n",alpha);
	printf("heat conductivity   (W/m K): %.3g \n",lambda);
	printf("Hello\n");
	return 0;

}
