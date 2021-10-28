// This is a rewrite of the original BRC program from G. Ahlers, but now it uses the CoolProp fluid
// library

// 		rho     =	density                  ( kg/m$3 )                   
//		alpha 	=	thermal expansion coeff	 (1/K )                       
//		comp	=	compressibility	 		 ( m^2/N )
//						(this is sometimes set to zero)               
//		lambda	= 	thermal conductivity 	 ( W/ m K )                    
//		kappa 	=	thermal diffusivity 	 ( m2/s )                         
//		nu 	    =	 shear viscosity 		 ( kg /s m )                      
//		cp 	    =   spec. heat at const. press. (J/kg K )
//
//	For the mixtures, the syntax is slightly different, e.g.:
// 
//	h2_xe(temp, press, x, &rho, &alpha, &comp, &lambda, &kappa, &nu, &cp, &psi, &lewis);
//
//	with
//		x	=	molar concentration of heavier component 
//		psi	=	separation ratio
//		lewis	=	Lewis number
//

#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include "libFluidPropC.h"
#include "header.def"
#include "header.h"

#define DIAMETER 0.1  		/* diameter of cylindrical cell in m */
#define HEIGHT 0.100 //9.964
#define CHOICE 2			/* 0 = fix dtc (onset), 1 fix height (onset), 2 fix height (turbulence) */

// Some material parameters for plate corrections
// TODO update units below
#define RHO_CU 8.92e3	/* g/cm^3 */
#define CP_CU 0.385e3	/* J/g/K */
#define RHO_AL 2.70e3	/* g/cm^3 */
#define CP_AL 0.9e3	/* J/g/K */

int main(){
    //double res{0.0};
	double T,P,temp, rho, alpha, comp, lambda, kappa, nu, cp;//, tscale;
	double temp1, rho1, alpha1, lambda1, kappa1, nu1, cp1;
	double temp2, rho2, alpha2, lambda2, kappa2, nu2, cp2;
	double gamma[5], corr[5][5];
	double p0, p1, p2, p3, p4;
	double q, rh20, rr20, l2;
	double epsa, epsb, epsr;
	double X, Y, Z, D0, a, b, epsT;
	double omeg, taylor, tv, U, rossby;
	double dtc, qcrit, sigma, sigma1, sigma2, dsigma, a_r, n_index, n_1, n_2, beta;
	double height, aspect, press;
	double g3, g3_tilde, eps_crit, qc, xi0, tau0, F_th;
	int gas, i, j;
	double psi, lewis, rayleigh, plate_ratio_CU, plate_ratio_AL;
	double x=1;
	double dtemp, volume, nusselt, lambda_eff, current, t_b, x_WL, alf_dt, drc, rc_gamma, a_g;
	char fluid[200],flag[100];
	
	sprintf(flag,"HEOS");

	if(CHOICE == 0)
		printf("\nDelta T_c is user specified. \nOutput for pattern formation near onset.\n");
	if(CHOICE == 1)
		printf("\nCell thickness is fixed by #define HEIGHT ... . \nOutput for pattern formation near onset.\n");
	if(CHOICE == 2)
		printf("\nCell thickness is fixed by #define HEIGHT ... . \nFor turbulent flows at large Rayleigh numbers.\n");
	printf("\nIf you would like a different option, change #define CHOICE \nin RBC.c and re-compile.\n\n");


	printf("Which fluid ? \n\t 0 = Air \t 1 = H2 \t 2 = He \t 3 = N2 \t 4 = CO2 \n\t 5 = Xe  \t 6 = SF6  \t 7 = SF6_crit\t 8 = Ethane\t 9 = Ethane_Crit ");
	printf("\n\t10 = 5CB \t11 = H2O \t12 = ACETONE \t13 = METHANOL ");
	printf("\t14 = ETHANOL \n\t15 = 2-PROPANOL\t16 = TOLUENE\t17 = GLYCEROL\t18 = triethyleneglycol ");
	printf("\n\t20 = H2-Xe \t21 = He-SF6 \t22 = Glycerol-Water\n");
	scanf("%d", &gas);

    if((gas < 10) || (gas > 20)){
            printf("mean temperature (deg C) and pressure (bar) ?   ");
            scanf("%lf %lf", &temp, &press);
			T = temp + T0;
			P  = press *1e5;
    }
    else{
            printf("mean temperature (deg C) ?   ");
            scanf("%lf", &temp);
			T = temp + T0;
            P = 1.0e5;
    }

	x = 0.; psi = 0.; lewis = 0.; 	/* for mixtures only : */
	if(CHOICE == 2){
		printf("Temperature difference ?  ");
		scanf("%lf", &dtemp);
	}

/* print out a general header */

	fn_header(gas, temp, press, x);

/* get properties at the mean temperature. NOTE: nu is SHEAR viscosity. */
	
	switch (gas){
		case 0:
			sprintf(fluid,"Air");
			break;
		case 1:
			sprintf(fluid,"Hydrogen");
			break;
		case 2:
			sprintf(fluid,"Helium");
			break;
		case 3:
			sprintf(fluid,"Nitrogen");
			break;
		case 4:
			sprintf(fluid,"CO2");
			break;
		case 5:
			sprintf(fluid,"Xenon");
			break;
		case 6:
			sprintf(fluid,"SF6");
			break;
		case 7:
			printf("SF6_crit  not yet implemented!!\n");
			break;
		case 8:
			sprintf(fluid,"Ethane");
			break;
		case 9:
			printf("C2H6_crit  not yet implemented!!\n");
			break;
		case 10:
			printf("5cb  not yet implemented!!\n");
			break;
		case 11:
			sprintf(fluid,"Water");
			break;
		case 12:
			sprintf(fluid,"Acetone");
			break;
		case 13:
			sprintf(fluid,"Methanol");
			break;
		case 14:
			sprintf(fluid,"Ethanol");
			break;
		case 15:
			printf("Isopropanol  not yet implemented!!\n");
			break;
		case 16:
			printf("Toluene  not yet implemented!!\n");
			break;
		case 17:
			printf("Glycerol  not yet implemented!!\n");
			break;
		case 18:
			printf("Triethyleneglycol  not yet implemented!!\n");
			break;
		case 19:
			printf("FC72 not yet implemented!!\n");
			break;
		case 20:
			printf("H2-Xe  not yet implemented!!\n");
			//printf("Molar concentration of the heavier component ?  ");
			//scanf("%lf", &x);
			break;
		case 21:
			printf("H2-SF6  not yet implemented!!\n");
			//printf("Molar concentration of the heavier component ?  ");
			//scanf("%lf", &x);
			break;
		case 22:
			printf("Weight concentration of Glycerol?  ");
			scanf("%lf", &x);
			sprintf(fluid,"MGL[%.2f]",x);
			sprintf(flag,"INCOMP");
			break;
		default:
			printf("no such fluid available\n");
			exit(0);

	}
	printf("T: %lf\tP: %lf\n",T,P);
	//shared_ptr<AbstractState> Water(AbstractState::factory("HEOS",fluid));
	getCoolProp(fluid, T, P, &rho, &alpha, &comp, &lambda, &kappa, &nu, &cp, &psi, &lewis,flag);
	//kappa  = lambda/(cp*rho);

	sigma = nu/kappa;				/* Prandtl number */
	plate_ratio_CU = sqrt(RHO_CU*CP_CU/(rho*cp*1.e-7));
	plate_ratio_AL = sqrt(RHO_AL*CP_AL/(rho*cp*1.e-7));

    a_g = alpha*(temp + T0)/(rho*cp);
    a_g = a_g*rho*G;

/* determine height and Delta T_c */

	if(CHOICE == 0){
		printf("Critical temperature difference (deg C) ?   ");
		scanf("%lf", &dtc);
		height = RC*kappa*nu/(rho*alpha*G*dtc);	
		height = pow(height, 1./3.);
	}
	else{
		height = HEIGHT;
		dtc = RC*kappa*nu/(rho*alpha*G*height*height*height);
printf("DTc = %.4e\n\n", dtc);
	}
	aspect = DIAMETER/height;
	//tscale = dtc/RC;

/* compute properties of turbulent state */

	if(CHOICE == 2){
		volume = 1.e-3*height*DIAMETER*DIAMETER*PI/4.;   /* cell volume in liter */
		rayleigh = G*dtemp*alpha*height*height*height/(kappa*nu);
/*
		nusselt = 0.3265*pow(rayleigh, 0.25)*pow(sigma, -1./12.);		
		nusselt += 0.002352*pow(rayleigh, 3./7.)*pow(sigma, -1./7.);
		nusselt = nusselt/(1. + 0.181*log10(aspect));*/

//		nusselt = NofR_Oregon(sigma, rayleigh);
//		nusselt = 0.137*pow(rayleigh, 0.300);
		nusselt = 0.104*pow(rayleigh, 0.312);
//		nusselt = 0.112*pow(rayleigh, 0.312);	
//      fn_nus(sigma, rayleigh, &nusselt, &reynolds);
		lambda_eff = lambda*nusselt;
	}

/* parameters needed near onset of RBC */

	qc = 3.117;
	xi0 = sqrt(0.148);
	tau0 = 1./(19.65*sigma/(sigma+0.5117));
		
/* thermal noise parameters, see Hohenberg and Swift, Phys. Rev. A 46, 4773 (1992). */

	F_th = (KB*(temp + T0)*rho/(height*nu*nu))*2.*sigma*qc/(xi0*tau0*RC);  /* Eq. 2.17 */
	g3 = 1.4;
	g3_tilde = (2./3.)*g3;
	eps_crit = 3.*g3_tilde*F_th/4.;
	eps_crit = pow(eps_crit, 2./3.);

/* get properties at hot and cold end */

	if((CHOICE == 0) || (CHOICE == 1)){
		temp1 = T + dtc/2.;
		temp2 = T - dtc/2.;
	}
	else if(CHOICE == 2){
		temp1 = T + dtemp/2.;
		temp2 = T - dtemp/2.;
	}
	else{
		printf("\nno such option.\n");
		exit(0);
	}

	getCoolProp(fluid, temp1, P, &rho1, &alpha1, &comp, &lambda1, &kappa1, &nu1, &cp1, &psi, &lewis,flag);
	getCoolProp(fluid, temp2, P, &rho2, &alpha2, &comp, &lambda2, &kappa2, &nu2, &cp2, &psi, &lewis,flag);
	
	//kappa1 = lambda1/(rho1*cp1);
	//kappa2 = lambda2/(rho2*cp2);
	sigma1 = nu1/kappa1;
	sigma2 = nu2/kappa2;

/* some refractive index stuff. Needed for shadowgraph sensitivity. */

    /* for refractive index calculation. OK only for CO2 and SF6.  */
        if(gas == 4)	
            a_r = 0.152;		/* CO2 */
        else if(gas == 6)
            a_r = 0.0777;		/* SF6 */
        else
            a_r = 0.;			/* all else */
	n_index = sqrt((1.+2.*a_r*rho)/(1.-a_r*rho));    /* refractive index */
	n_1 = sqrt((1.+2.*a_r*rho1)/(1.-a_r*rho1));
	n_2 = sqrt((1.+2.*a_r*rho2)/(1.-a_r*rho2));
	beta = (n_1 - n_2)/dtc;


/* non-Boussinesq things */
	
	if(CHOICE == 0 || CHOICE == 1)
		dsigma = 2./(sigma1+sigma2)*(sigma2-sigma1)/dtc;
	else if(CHOICE == 2)
		dsigma = 2./(sigma1+sigma2)*(sigma2-sigma1)/dtemp;
	
	gamma[0] = -(rho1 - rho2)/rho;
	gamma[1] = (rho1*alpha1-rho2*alpha2)/(2.*rho*alpha);
	gamma[2] = (nu1/rho1 - nu2/rho2)/(nu/rho);
	gamma[3] = (lambda1 - lambda2)/lambda;
	gamma[4] = (cp1 - cp2)/cp;

	/* Busse coefficients. These are not used. We use the new Pesch coefficients.*/
/*	 
	printf("NOB coeffs from Busse 1967\n");
	p0 = 2.676 - 0.1258/sigma;
	p1 = -6.603 - 0.5023/sigma;
	p2 = 2.755;
	p3 = 2.917 - 0.5023/sigma;
	p4 = -6.229 + 0.2512/sigma;
*/
	/* Pesch coefficients */
	printf("NOB coeffs from BPA00, Annu. Rev. Fluid Mech. 32, 709 (2000), Sect. 6.5\n");
	p0 = 2.676 - 0.361/sigma;
	p1 = -6.631 - 0.772/sigma;
	p2 = 2.765;
	p3 = 9.540;
	p4 = -6.225 + 0.3857/sigma;

	l2 = 0.29127 + 0.08147/sigma + 0.08933/(sigma*sigma);
	rh20 = 0.89360 + 0.04959/sigma + 0.06787/(sigma*sigma);
	rr20 = 0.69942 - 0.00472/sigma + 0.00832/(sigma*sigma);
	q = gamma[0]*p0 + gamma[1]*p1 + gamma[2]*p2 + gamma[3]*p3 + gamma[4]*p4;

/* these formulas are from Busse 1967*/
	epsa = -q*q/(4.*rh20*RC); 
	epsb = q*q*(9.*rh20 - 3.*l2)/(l2*l2*RC);		
	epsr = q*q*3.*rr20/(l2*l2*RC);

/* the following is from BdAC91 [PRL 67, 3078 (1991)] */
	X = 1.39892 - 0.00944/sigma + 0.01665/(sigma*sigma);
	Y = 1.98148 + 0.15349/sigma + 0.19529/(sigma*sigma);
	a = q*sqrt(6./(RC*X));
	b = Y/X;
	D0 = 1. + b - 2.*b*b;
	Z = a*b/D0 - (a/fabs(a))*sqrt(a*a*b*b/(D0*D0) + a*a/(2.*D0));
		/* Z is corrected ala footnote on p. 744 of BPA00 */ 
	epsT = a*Z+(2.*b+1)*Z*Z; 	

/* the following formulas are from BdAC91 [PRL 67, 3078 (1991)], and give the same results */
/*  as the Busse formulas used above for epsa, etc. 
	epsa2 = -a2/(4.+8.*b);
	epsr2 = a2/((1.-b)*(1.-b));
	epsb2 = a2*(2.+b)/((1.-b)*(1.-b));
*/

/* second order non-OB shift of R_c */
/* coefficients corr[i][j] from W. Pesch, private commun. */

        //bzero(corr,25*sizeof(double));
        corr[0][0] = 10.122;
        corr[1][0] = -12.504;
        corr[1][1] = -12.984;
        corr[2][0] = 0.;
        corr[2][1] = 76.833;
        corr[2][2] = -88.059;
        corr[3][0] = -0.00157;
        corr[3][1] = 153.738;
        corr[3][2] = -76.844;
        corr[3][3] = -76.923;
        corr[4][0] = 6.252;
        corr[4][1] = -51.861;
        corr[4][2] = 38.417;
        corr[4][3] = 32.436;
        corr[4][4] = -3.246;
        
        drc = 0.;
        for(i = 0; i <= 4; i++){
            for(j = 0; j <= i; j++){
                drc += gamma[i]*corr[i][j]*gamma[j];
            }
        }
        rc_gamma = 1707.8 + drc;

/* correct dT_c for NOB effects */
/* Note that the gammas's are based on R_c = 1707.8 when dT_c is calculated. */

	dtc = dtc*rc_gamma/1707.8;

/* Boundary-layer calculations, Eq. 8 of Wu and L., Phys. Rev. A 43, 2833 (1991) */

	if(CHOICE == 2){
	/* I will use an estimate of the average temperature of each BL to get the properties : */
		temp1 = temp + dtemp/4.;
		temp2 = temp - dtemp/4.;
		x_WL = alpha1*(nu2/rho2)*kappa2/(alpha2*(nu1/rho1)*kappa1);
		x_WL = pow(x_WL, 1./3.)*lambda1/lambda2;
		alf_dt = alpha*dtemp;		/* NOB criterion used by the Oregon group */
	}

/* End of non-Boussinesq */

// Addition by S. Weiss (January 2007)
/* Lets calculate the real temperature difference between the top of top plate and bottom of bottom plate */

/* common to all options: cell geometry and fluid props at mean temp. */

	printf("\nHeight = %.4f m\t\tAspect Ratio = %.3f\n\n", height, aspect);

	printf("density (kg/m^3)                 = %.5g \n",rho);	
	printf("isothermal expansion coeff (1/K) = %.6g \n",alpha);
	printf("comp (m^2/N)                     = %.5g \n",comp);
	printf("heat conductivity (W/K m)        = %.5g \n",lambda);
	printf("thermal diffusivity (m^2/s)      = %.5g \n",kappa);
	printf("kinematic viscosity (m^2/s)      = %.5g \n",nu);
	printf("isobaric heat capacity (J/K kg)  = %.5g \n",cp);
	printf("Prandtl number                   = %.3f\n\n", sigma);
	printf("\n");	
	printf("dynamic viscosity (Pa s)         = %.5g \n",nu*rho);
	printf("[(rho*Cp)_Cu/(rho*Cp)]^(1/2) = %.3f\n", plate_ratio_CU);
	printf("[(rho*Cp)_Al/(rho*Cp)]^(1/2) = %.3f\n", plate_ratio_AL);
	if(gas >= 20)
		printf("x = %.3f\t\tPsi = %.3f\t\t\tLewis = %.3f\n\n", x, psi, lewis);
	printf("\n");

/* Output relevant to RBC onset */

	if(CHOICE == 0 || CHOICE == 1){
		qcrit = dtc*lambda*DIAMETER*DIAMETER/4.*PI/height*1.e-7; 
		tv = height*height*rho/nu;
		omeg = 2.*3.14159*tv*0.25;
		taylor = 4.*omeg*omeg;
                printf("gamma_i :	");
                for(i = 0; i < 5; i++)
                    printf("%.4f   ", gamma[i]);
                printf("\n");
		printf("F_th = %.3e\t\teps_crit = %.3e\n", F_th, eps_crit);
		printf("Q_Busse = %.3f\t\tEpsa  = %.3e\t\tEpsr  = %.3e\n", q, epsa, epsr);
		printf("EpsT' = %.3e\t\tEpsb  = %.3e\n", epsT, epsb);
                printf("R_c^gamma = %.1f	dR_c = %.3f	a_g = %.3e K/cm\n", rc_gamma, drc, a_g);
		printf("delta tc = %.4e C		Qcrit = %.4e W\n",dtc, qcrit);

		printf("visc. tv = %.3f s		therm. tv = %.3f s\n\n",tv, tv*sigma);	  
		printf("1/4 Hz yields Omega = %.3e and\tTaylor = %.3e\n", omeg, taylor);
		if((gas == 4) || (gas == 6)){
//			sens = -9.44e4*beta*height*dtc;
			printf("Refractive index n = %.5f   dn/dT = %.3e\n\n", n_index, beta);
		}

	}

/* Output relevant to turbulent RBC */

	else if(CHOICE == 2){
		printf("Q_Busse = %.3f\t\tx_WL = %.4f\t\talpha*DeltaT = %.4f\n", q, x_WL, alf_dt);
		qcrit = dtc*lambda*DIAMETER*DIAMETER/4.*PI/height*1.e-7; 
		printf("delta tc = %.4e C		Qcrit = %.4e W\n",dtc, qcrit);
		current = dtemp*lambda*nusselt*height*aspect*height*aspect/4.*PI/height*1.e-7; 
		printf("Rayleigh = %.2e\t\tNusselt = %.2f\nPrandtl = %.3f\t\t\n", rayleigh, nusselt, sigma);
		printf("Q = %.3f W\t\t\tVolume = %.3f liter\n", current, volume);

		tv = height*height*rho/nu;
		t_b = height*height*rho*cp/lambda_eff;

		printf("visc. tv = %.3f s\ttherm. tv = %.3f s\tlambda_eff = %.3e ergs/s cm K\n\n",tv, t_b, lambda_eff);
		
//		double q_over_dts = 33.3397 + 1.0762e-2*(temp + dtemp/2.);
		double q_over_dts = 33.3397 + 1.0762e-2*(temp + dtemp/2.);
		double dts = current/q_over_dts;
		printf("Ttop = %.3f\tT_bot_t = %.3f\tT_bot_b = %.3f\n\n", temp-dtemp/2., temp+dtemp/2., temp+dtemp/2.+dts);
		
		tv = height*height*rho/nu;
		omeg = 2.*3.14159*tv*0.25;
		taylor = 4.*omeg*omeg;
		U = sqrt(G*alpha*dtemp*height);	// "free-fall velocity" in cm/sec
		rossby = U/(2.*0.5*PI*height);			// for Omega = 1 Hz
		printf("1/4 Hz yields Omega = %.3e and\tTaylor = %.3e\n", omeg, taylor);
		printf("1/4 Hz yields Rossby = %.4f and\tU = %.3f\n", rossby, U);
	}

	else{
		printf("\nno such option.\n");
		exit(0);
	}

	return 0;


}


void fn_header(int gas, double temp, double press, double x)
{

	printf("gas = %d\n", gas);

	switch (gas){
		case 0:
			printf("\n\nAir\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 1:
			printf("\n\nH2\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 2:
			printf("\n\nHe\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 3:
			printf("\n\nN2\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 4:
			printf("\n\nCO2	\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 5:
			printf("\n\nXe\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 6:
			printf("\n\nSF6\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 7:
			printf("\n\nSF6 near the critical point\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 8:
			printf("\n\nEthane\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 9:
			printf("\n\nEthane near the critical point\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\n\n",temp, press);
			break;
		case 10:
			printf("\n\n5CB\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 11:
			printf("\n\nH2O	\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 12:
			printf("\n\nACETONE\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 13:
			printf("\n\nMETHANOL\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 14:
			printf("\n\nETHANOL\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 15:
			printf("\n\n2-PROPANOL\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 16:
			printf("\n\nTOLUENE\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 17:
			printf("\n\nGLYCEROL\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 18:
			printf("\n\nTRIETHYLENEGLYCOL\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 19:
			printf("\n\nFC72\t\tT_bar = %6.3f deg C\n\n",temp);
			break;
		case 20:
			printf("\n\nH2-Xe\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\t\tX_Xe = %6.3f\n\n",temp, press, x);
			break;
		case 21:
			printf("\n\nHe-SF6\t\tT_bar = %6.3f deg C\t\tP = %6.2f bar\t\tX_SF6 = %6.3f\n\n",temp, press, x);
			break;
		case 22:
			printf("\n\nGlycerol-Water\t\tT_bar = %6.3f deg C\tGlycerol = %6.3f\n\n",temp, x);
			break;
		default:
			printf("no such fluid\n");
			exit(0);
		}
}


double NofR_Oregon(double sig, double r)
{
	/* A fit to the Oregon data, after correction for wall conductance using Model 2. */
	double n, ln, lr, a, b;
	a  = -1.140182;
	b  = 0.324721;
	lr = log10(r);
	ln = a + b*lr;
	n  = pow(10., ln);
	return(n);
}
