#include <stdio.h>
#include "libFluidPropC.h"
#include "CoolPropLib.h"
#include "AbstractState.h"
#include "crossplatform_shared_ptr.h"

extern "C" {
	int getCoolProp(char * fluid, double T, double P, double x,double *prho, double *palpha, double  *pcomp, double *plambda, double *pkappa, double *pnu, double  *pcp, double *ppsi, double *plewis ){
		//*prho    = PropsSI("D", "T", T, "P", P, fluid);
		//*plambda = PropsSI("L", "T", T, "P", P, fluid);
		//*pcp     = PropsSI("CPMASS", "T", T, "P", P, fluid);
		//*palpha  = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T, "P", P, fluid);
		//*pnu     = PropsSI("VISCOSITY", "T", T, "P", P, fluid)/(*prho);
		//*pkappa  = *plambda/((*pcp)*(*prho));
		
		using namespace CoolProp;
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
		printf("P: %lf\tT: %lf\trho: %lf\n",P,T,*prho);
		return 0;
	}
}
