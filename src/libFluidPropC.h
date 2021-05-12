#ifndef __LIBCOOLPROPC_H
#define __LIBCOOLPROPC_H

#ifdef __cplusplus
extern "C" {
#endif
	int getCoolProp(char * fluid, double T, double P, double *prho, double *palpha, double  *pcomp, double *plambda, double *pkappa, double *pnu, double  *pcp , double *ppsi, double *plewis);

#ifdef __cplusplus
}
#endif
#endif
