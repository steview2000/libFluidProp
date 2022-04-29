#include <stdio.h>
#include "libFluidPropC.h"
#include "CoolPropLib.h"
#include "AbstractState.h"
#include "crossplatform_shared_ptr.h"

extern "C" {
	int getCoolProp(char * fluid, double T, double P, double *prho, double *palpha, double  *pcomp, double *plambda, double *pkappa, double *pnu, double  *pcp, double *ppsi, double *plewis,char *flag ){
		
		using namespace CoolProp;
		double Tcrit,Pcrit,Pvap;
		//int phase;
		//double rho,rho1,rho2,T1,T2;
		//char fullfluid[100];
		
		shared_ptr<AbstractState> fluid_PTR(AbstractState::factory("HEOS",fluid));	
		if (strcmp(flag,"SLOW")==0){
			// get critical point
			Tcrit = fluid_PTR->keyed_output(iT_critical)	;
			Pcrit = fluid_PTR->keyed_output(iP_critical)	;
			printf("Tcrit: %.4lfK,Pcrit: %.4lfpa\n",Tcrit,Pcrit);
			if (T<Tcrit){ 
				if (P<Pcrit){
					printf("Fluid-string: %s %s\n",flag,fluid);
					Pvap = PropsSI("P", "T", T, "Q", 0, fluid);
					printf("Pvap: %lf\n",Pvap);
					if (Pvap>P){	
						printf("Phase: Gas\n");
						fluid_PTR->specify_phase(iphase_gas)	;
					}
					else printf("Phase: Liquid\n");
				};
			}
			else printf("Supercritical\n");
			//if (strcmp(flag,"GAS")== 0) fluid_PTR->specify_phase(iphase_gas);
			//else if (strcmp(flag,"LIQUID")== 0) fluid_PTR->specify_phase(iphase_liquid);
			//else if (strcmp(flag,"SUPERCRITICAL_GAS")== 0) fluid_PTR->specify_phase(iphase_supercritical_gas);
			//else if (strcmp(flag,"SUPERCRITICAL_LIQUID")== 0) fluid_PTR->specify_phase(iphase_supercritical_liquid);
		}

		fluid_PTR->update(PT_INPUTS, P, T);
		*prho    = fluid_PTR->rhomass()	;
		*pcp     = fluid_PTR->cpmass()	;
		*palpha  = fluid_PTR->isobaric_expansion_coefficient();
		*plambda = fluid_PTR->conductivity();
		*pnu     = fluid_PTR->viscosity()/(*prho);
		*pcomp   = fluid_PTR->isothermal_compressibility();
		*pkappa   = *plambda/(*pcp * *prho);

		// TODO add the pcomp and ppsi and plewsi
		*ppsi   = 0;
		*plewis = 0;
		//printf("P: %lf\tT: %lf\trho: %lf\n",P,T,*prho);
		//sprintf(fullfluid,"HEOS::%s",fluid);
		//printf("SLowFLUID\n");
		//T1 = T*1.0001;
		//T2 = T*0.9999;

		//rho1    = PropsSI("D", "T", T1, "P", P, fullfluid);
		//rho2    = PropsSI("D", "T", T2, "P", P, fullfluid);
		//rho = 0.5*(rho1+rho2);
		////printf("P: %lf\tT1: %lf\t T2: %lf\n",P,T1,T2);
		////printf("rho1: %lf\t rho2: %lf\n",rho1,rho2);
		//*prho     = rho;
		//*plambda = PropsSI("L", "T", T, "P", P, fullfluid);
		//*pcp     = PropsSI("CPMASS", "T", T, "P", P, fullfluid);
		//*palpha  = (rho1-rho2)/(rho*(T2-T1));//PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T, "P", P, fluid);
		//*pnu     = PropsSI("VISCOSITY", "T", T, "P", P, fullfluid)/(rho);
		//*pkappa  = *plambda/((*pcp)*(*prho));
		//*pcomp   = 0;
		//printf("P: %lf\tT: %lf\tnu: %lf\n",P,T,*pnu);
		return 0;
	}
}
