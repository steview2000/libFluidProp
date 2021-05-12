#include <stdio.h>
#include "libFluidPropC.h"
#include "include/CoolPropLib.h"

extern "C" {
	int getCoolProp(char * fluid, double T, double P, double *prho, double *palpha, double  *pcomp, double *plambda, double *pkappa, double *pnu, double  *pcp, double *ppsi, double *plewis ){
		*prho    = PropsSI("D", "T", T, "P", P, fluid);
		*plambda = PropsSI("L", "T", T, "P", P, fluid);
		*pcp     = PropsSI("CPMASS", "T", T, "P", P, fluid);
		*palpha  = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", T, "P", P, fluid);
		*pnu     = PropsSI("VISCOSITY", "T", T, "P", P, fluid)/(*prho);
		*pkappa  = *plambda/((*pcp)*(*prho));

		// TODO add the pcomp and ppsi and plewsi
		*pcomp = 0;
		*ppsi = 0;
		*plewis = 0;
		printf("P: %lf\tT: %lf\trho: %lf\n",P,T,*prho);
		return 0;
	}
}
