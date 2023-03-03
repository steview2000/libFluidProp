
/* -------  general constants ---------- */

#define ZERO_C 	273.15			/* 0 deg. C in K */
#define T0		273.15			/* same as ZERO_C */
#define PI 		3.14159263
#define GAMA	0.577215
#define R 		8314.33    		/*  J/kmole*K  */
#define NA		6.0221367e+23	/* Avogadro's number */
#define KB 		1.38027E-23		/* Boltzmann's constant erg/K */
#define K		1.38027e-23		/*   J/K   */
#define RG 		8.314462			/* gas constant in J/mole K */
#define RC		1707.8			/* Crit. Rayleigh number */
#define G 		9.81	/* local gravitational acceleration  in m/sec^2 */
#define PCONV	14.504			/* psi per bar */

/* ------ fluid number ----*/
#define AIR 0
#define HYDROGEN 1 
#define HE 2 
#define N2 3 
#define CO2 4 
#define XENON 5 
#define SF6 6

/*-------- refractive-index virial coefficients -------*/

#define AR_HE	5.2e-7      /* unit:  m^3/mol   */
#define AR_CO2	6.65e-6
#define AR_SF6	1.134e-5
#define BR_HE	-6.0e-14    /* unit:  m^6/mol^2   */
#define BR_CO2	3.0e-12
#define BR_SF6	2.8e-11
#define AR_XE   1.036e-5 /* at 298 K */	
#define BR_XE	2.58e-11 /* at 298 K deg C */
#define AR_H2	2.0713e-6      /* unit:  m^3/mol   */
#define BR_H2	1.3e-13      /* unit:  m^6/mol^2   */

/* -------  molecular mass ----------- */

#define MW_H2 			2.0156   		/*  g/mol    */
#define MW_HE 			4.0026
#define MW_N2 			28.0134 
#define MW_CO2 			44.01
#define MW_XE 			131.29
#define MW_SF6			146.066 
#define MW_ACETONE		58.08
#define MW_d_ACETONE	64.13		/* deuterated Acetone */
#define MW_METHANOL	32.042
#define MW_TOLUENE		92.141 
#define MW_ETHANOL		46.069
#define MW_PROP			60.096
#define MW_GLYCEROL	92.095
#define MW_5CB 			249.4
#define MW_TRIETHYLENEGLYCOL 150.18

#include <stdio.h>
#include <stdlib.h>                                         
#include <math.h>

void fn_header(int, double, double, double);
void get_prop(int, double, double, double, double *, double *, double *, double *, double *,  double *, double *, double *, double *);
void cdo(double, double, double *, double *, double *, double *, double *, double *, double *);
void fn_5cb(double, double, double *, double *, double *, double *, double *, double *, double *);
void water2(double, double, double *, double *, double *, double *, double *, double *, double *);
void acetone(double, double, double *, double *, double *, double *, double *, double *, double *);
void nitrogen(double, double, double *, double *, double *, double *, double *, double *, double *);
void helium(double, double, double *, double *, double *, double *, double *, double *, double *);
void sf6(double, double, double *, double *, double *, double *, double *, double *, double *);
void xenon(double, double, double *, double *, double *, double *, double *, double *, double *);
void ethanol(double, double, double *, double *, double *, double *, double *, double *, double *);
void methanol(double, double, double *, double *, double *, double *, double *, double *, double *);
void isopropanol(double, double, double *, double *, double *, double *, double *, double *, double *);
void toluene2(double, double, double *, double *, double *, double *, double *, double *, double *);
void glycerol(double, double, double *, double *, double *, double *, double *, double *, double *);
void sf6_crit(double, double, double *, double *, double *, double *, double *, double *, double *);
void h2_xe(double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *);
void he_sf6(double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *);
double NofR_GL(double, double);
double NofR_Oregon(double, double);
double refIndexTayag(double ,double ,int );
double calc_f(double x);
double calc_g(double x);
double grossmannLohse(double,double);
int printUsage();
