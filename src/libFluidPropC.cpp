#include <stdio.h>
#include "libFluidPropC.h"
#include "CoolPropLib.h"
#include "AbstractState.h"
#include "crossplatform_shared_ptr.h"

extern "C" {
	int getCoolProp(char * fluid, double T, double P, double *prho, double *palpha, double  *pcomp, double *plambda, double *pkappa, double *pnu, double  *pcp, double *ppsi, double *plewis,char *flag ){
		
		using namespace CoolProp;
		double rho,rho1,rho2,T1,T2;
		char fullfluid[100];


		//printf("Fluid-string: %s %s\n",flag,fluid);

		if (strcmp(flag,"INCOMP")== 0){
			sprintf(fullfluid,"INCOMP::%s",fluid);
			//printf("INCOMPRESSIBLE FLUID\n");
			T1 = T*1.01;
			T2 = T*0.99;

			rho1    = PropsSI("D", "T", T1, "P", P, fullfluid);
			rho2    = PropsSI("D", "T", T2, "P", P, fullfluid);
			rho = 0.5*(rho1+rho2);
			//printf("T1: %lf\t T2: %lf\n",T1,T2);
			//printf("rho1: %lf\t rho2: %lf\n",rho1,rho2);
			*prho     = rho;
			*plambda = PropsSI("L", "T", T, "P", P, fullfluid);
			*pcp     = PropsSI("CPMASS", "T", T, "P", P, fullfluid);
			*palpha  = (rho1-rho2)/(rho*(T2-T1));//PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T, "P", P, fluid);
			*pnu     = PropsSI("VISCOSITY", "T", T, "P", P, fullfluid)/(rho);
			*pkappa  = *plambda/((*pcp)*(*prho));
			*pcomp   = 0;
			//printf("P: %lf\tT: %lf\tnu: %lf\n",P,T,*pnu);
		}
		else{
			shared_ptr<AbstractState> fluid_PTR(AbstractState::factory("HEOS",fluid));	

			fluid_PTR->update(PT_INPUTS, P, T);
			*prho    = fluid_PTR->rhomass()	;
			*pcp     = fluid_PTR->cpmass()	;
			*palpha  = fluid_PTR->isobaric_expansion_coefficient();
			*plambda = fluid_PTR->conductivity();
			*pnu     = fluid_PTR->viscosity()/(*prho);
			*pcomp   = fluid_PTR->isothermal_compressibility();

			// TODO add the pcomp and ppsi and plewsi
			*ppsi   = 0;
			*plewis = 0;
			//printf("P: %lf\tT: %lf\trho: %lf\n",P,T,*prho);
		}
		return 0;
	}
}
