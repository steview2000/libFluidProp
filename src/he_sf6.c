/*****************************************************************************
hsf.c

Mixture of helium and sulfur hexafluoride.  SI unit.

The program is valid for T > 208 K and quite low density (ro << B/C).  The 
zero-density formulas are taken from J. Bzowski et al, J. Phys. Chem. Ref. 
Data, V.19, 1179 (1990), which is based on the kinetic theory of corresponding 
states.  The parameters used in the equations had been carefully calculated by 
the authors from experimental data available to them.

The interaction virial coefficients B_12 (nor B_mix) are larger than the 
measured ones by about 0.02 m^3/kmol, so that I corrected this 'error' by 
substracting 0.018 m^3/kmol from calculated B_12. 

11/14/94, by Jun Liu.

(1)  Add the binary diffusion coefficient D and thermal diffusion factor 
 k_t to the program.  The formulas are also from Bzowski et al.'s paper (but 
 eq. C5 in the paper lost a factor 6.0 in front of Cstar_ij).  For k_t, A few 
 available experimental points are about 10% to 15% larger then the results of 
 Bzowski's program, so I also correct this 'error'. 

 Both D and k_t are also presure dependent [T. N. Bell et al, Chem. Phys. Lett. 
 V.45, 445 (1977) and R.D. Trengove et al. Physica A V.108, 502 (1981)]. The 
 experimental measurements only went up to 9 atm for D and 5 atm for k_t.  
 However, I use the second virial coefficients from these experiments for our 
 calculations.  D changes a couple of percent from above calculations over
 10 atm but k_t can change quite a lot.  

(2) Give a very rough estimation of dynamic viscosity as a function of 
 pressure.  Since the dynamic viscosity of helium does not change with the
 pressure, the variation of the mixture viscosity is mainly contributed by SF6, 
 so that we copy the density-dependent part of the viscosity of SF6 into this 
 program.  Such an estimation is OK considering that the viscosity of SF6 does 
 not vary very much with pressure (9.2% from 1 to 18 bar for pure SF6).

3/31/95, by Jun Liu 

The heat capacity is estimated by the summation of individual components.
Since the heat capacity of SF6 is pressure-depnedent, we need to figure out 
its partial pressure when we calculate the contribution of SF6. There are no 
experimental data on heat capacity of the mixture.

To calculate the Dufour number, we need to know the chemical potential. Here I  
use the formula for the ideal gas to give an rough estimation.  Considering we 
usually run our experiments at low concentration of SF6,and Helium is quite 
close to ideal, such an estimation should not be bad.  

4/10/95, by Jun Liu

  One important point should be clear that: in the governing equations of a
  binary gas mixture (Landau & Lifshitz, Fluid Mechanics; Hort et al. Phys.
  Rev. A45, 3737 (1992)), the concentration is MASS concentration.  Therefore
  k_t in these equations is related to mass concentration instead of the mole
  fraction usually used by people who measure the thermal diffusion factor. 

4/23/95, by Jun Liu

 Include higher (C and D) virial coefficients which are based on those of pure gases.

10/10/95, by Jun Liu

  Add refractive index. For more information see gasindex.c

10/13/95, by Jun Liu

Thermal conductivity is now based on a model proposed for polyatomic gases by 
J. Kestin and W. A. Wakeham, Int. J. Thermophys. V.4, 295 (1983). They employed 
measured physical properties whenever possible.  A. A. Clifford, J. Kestin and 
W. A.  Wakeham, Ber. Bunsenges. Phys. Chem. V.85, 385 (1981) also showed the method
to be better than direct calculations from the kinetic theory.  However, the
astar and bastar are still from the kinetic theory.

It turns out that the thermal conductivity given by this method is closer to my
measurement than the value given by a method due to Cheung et al. AIChE J. V.8, 221(1962).  

9/25/96, by Jun Liu

Now we calculate the pressure dependence of binary diffusion constant D using an extended model of the Thorne-Enskog theory proposed in T.N. Bell, I. R. Shankland and P. J. Dunlop, Chem. Phys. Lett. V.45, 445 (1977).

10/17/96, Jun Liu


  Limitation:  0 C < T < 100 C (for large T may be OK, but I did not check
				everything, so some results may be bad.)
   	       0 < C_sf6 <= 1  [When C_sf6 = 0 (pure helium), some quantities
				cannot be calculated, try no-zero but very 	
				small C_sf6.  For pure gas properties, you are 
				recommanded to use helium-p and sf6-p.]
               P:  If C_sf6 is quite large, P cannot be large, e.g. for
				C_sf6=1, P cannot be large then 20 bar.  
				Acctually pure SF6 is liquified when p is 
				above 22 bar.  I do not know exactly P's
				upper limit for this program.  I do not
				recommand this program for P > 60 bar.

  Precision (based on available data and not very extreme conditions!  The 
		errors are only rough estimation, not a strict statistics.  The 	
		errors become larger at higher pressure.):

	dynamic viscosity (pressure)     <  				2-3%
	second virial coefficient	 <	   			2-4%
        density 			 <	   			2-7%
	thermal conductivity		 for pure liquid     <		2-5%	
					 not known for mixture (guess < 10%?)
        diffusion			 <	 			2-3%
	thermal diffusion ratio		 a few experimental points      < 5%
	heat capacity			 for pure liquid     <		2-3%
					 nor known for mixture (guess < 10%?)

Last modified by Jun Liu on 10/25/96.
******************************************************************************/
/*  ----  Header -----  */

#include "header.h"

/* ------- parameters from Bzowski et al.'s paper --------- */

#define SIGMA_HSF6 	0.4298e-9 		/* unit: m */
#define SIGMA_SF6 	0.5252e-9			
#define SIGMA_HE	0.2610e-9
#define EK_HSF6		19.24			/* unit: Kelvin */
#define EK_SF6		207.7
#define	EK_HE		10.40
#define ROSTAR_HE	0.0797			/*  dimensionless  */
#define ROSTAR_SF6	0.050
#define ROSTAR_HSF6	0.0548
#define VSTAR_HE	8.50e+5			/*  dimensionless  */
#define VSTAR_SF6	4.067e+8
#define VSTAR_HSF6	7.590e+7

/* ------ gas refractive index ------------*/

#define AR_HE	5.2e-7      /* unit:  m^3/mol   */
#define AR_SF6	1.134e-5
#define BR_HE	-6.0e-14    /* unit:  m^6/mol^2   */
#define BR_SF6	2.8e-11

  /*  The above coefficients are from H.J. Achtermann et al. J. Chem. Phys.
	Vol.94, 5669 (1991), A.D. Buckingham et al. Proc. R. Soc. Lond. A. 
	Vol.336, 275 (1974) and J.M. St-arnaud et al. J. Chem. Phys. Vol.71, 
	4951 (1979).	*/

/* ------ controling number used in the program  ----- */

#define NP 		100
#define DDT		0.001        	/*  temperature increment  (K) */
#define DDP 		1.0		/*  pressure increment (Pa) */
#define DDC		1.0e-5		/*  molar concentration increment */
#define JMAX		26

he_sf6(tc, pp, x_sf, pro, palp, pcomp, plambmix, pkapap, petamix, pcpmix, ppsi, pL)
double tc, pp, x_sf, *pro, *palp, *pcomp, *plambmix, *pkapap, *petamix, *pcpmix, *ppsi, *pL;
{

double lambmix,lambda_tp,kapap;		/*  thermal conductivity and 
								diffusivity */
double cpmix, cpmixmass;   			/*  isobaric heat capacity, per mole or
								per unit mass.  unit: kJ/kmol*K
								or kJ/kg*K  */
double etamix, eta_p, eta_tp, nu;      	 /*  viscosity  */
double bmix, cmix, dmix;				/*  second virial coefficient */
double Dh, Dsf, Dmix, Dmp;			/*  diffusion constants */
double ktmix, ktmp, rktmix, rktp, rckt;  
								/*  thermal diffusion ratio, molar
 					    			concentration related, the 
					    			symbol with "r" means reduced, 
 					    			with "c" means corrected.  */
double ktmixmass, ktmpmass;		/*  mass concentration related */
double dmudcm;
double pbar;						/* pressure in bar */
double rrsf, rrdsf, rrhe;				/*  density of sf6 */
double tk;    						/*  temperature etc */
double ppa;	      					/*  pressure */
double ro;							/*  density  */
double alp;						/*  thermal expansion coeff.  */
double x_h;						/*  molar concentration */
double cmass_sf, dcdx, dxdc;			/*  mass concentration */
double Pr, L, Q, Rf;					/*  Prandtl, Lewis and Dufour number */
double beta, betamass, ksi;			/*  solutal, baric expansion coeff.  */
double psi, psip;					/*  separation ratio */	
double s_c, s_t, n, n_c_hat, n_t_hat;	
double deltaT, thick;

	//pbar = pp/PCONV;				/* convert pressure from psi to bar */ 
	pbar = pp;				/* not necessar anymore */ 

/* conversion */

	tk = tc + ZERO_C;
	ppa = pbar*1.0e+5;
	x_h = 1.0 - x_sf;

/* call subroutines to do calcualtion */

 	virial(x_sf, tk, &bmix, &cmix, &dmix);
	vislamb(x_sf, tk, ppa, bmix, &etamix, &lambmix, &cpmix);
	diff(x_sf, tk, ppa, &Dh, &Dsf, &Dmix, &rktmix);
	rhoalp(x_sf, tk, ppa, &ro, &alp, &beta, &ksi, &n, &n_c_hat, &n_t_hat);
	dktpress_he_sf6(ppa, x_sf, ro, tk, Dmix, rktmix, &Dmp, &ktmix, &ktmp, &rckt, &rktp);
	chemical_he_sf6(tk, x_sf, &dmudcm);	

	cmass_sf = x_sf*MW_SF6/(x_sf*MW_SF6 + x_h*MW_HE);
	rrsf = ro*cmass_sf;
	rrhe = ro - rrsf;

  /* estimate the variation of viscosity with pressure, which is contributed by SF6 (see program sf6-p.c). The contribution of helium is so weak that it is neglected. */


	rrdsf = rrsf/6.6159;					/*  reduced density  */
	eta_p = 1.3865e-2*rrdsf + 2.85099e-3*rrdsf*rrdsf - 4.16848e-5*pow(rrdsf, 3) + 6.446345e-7*pow(rrdsf, 4);
	eta_p = eta_p - 4.99704e-9*pow(rrdsf, 5) + 1.970048e-11*pow(rrdsf, 6) - 2.961807e-14*pow(rrdsf, 7);
	
	eta_p = eta_p*1.0e-6;       			/* get unit right */   
	eta_tp = etamix + eta_p;

	nu = eta_tp/ro;						/* kinematic viscosity */

  /* estimate the variation of thermal conductivity with pressure beased on the two pure components. */

    /* Helium part  */

	lambda_tp = lambmix + 2.33e-4*rrhe + 2.39e-8*rrhe*rrhe;
    
    /* SF6 part  */

	lambda_tp = lambda_tp + 1.1e-5*rrsf - 7.615e-9*rrsf*rrsf + 2.887e-10*rrsf*rrsf*rrsf;

  /* convert beta from the derivative to molar concentration to the 
	derivative to mass concentration   */

	dcdx = (x_sf*MW_SF6 + x_h*MW_HE)*(x_sf*MW_SF6 + x_h*MW_HE);
	dcdx = MW_SF6*MW_HE/dcdx;

        dxdc = ((1.0 - cmass_sf)*MW_SF6 + cmass_sf*MW_HE);
	dxdc = dxdc*((1.0 - cmass_sf)*MW_SF6 + cmass_sf*MW_HE);
	dxdc = MW_SF6*MW_HE/dxdc;

	betamass = beta*dxdc;

        cpmixmass = cpmix/(x_h*MW_HE + x_sf*MW_SF6);
						/* cpmixmass unit: kJ/kg*K  */

/* dimensionless numbers */

  /*  Prandtl number  */

	Pr = 1000.0*cpmix*eta_tp/(lambda_tp*(x_h*MW_HE + x_sf*MW_SF6));
			                       /* cpmix unit:  kJ/kmole*K   */

  /*   Lewis number   */

	kapap = lambda_tp*(x_h*MW_HE + x_sf*MW_SF6)/(1000.0*cpmix*ro);
        L = Dmp/kapap;

   /*  Separation ratio.  Since the ktmix calculated above are related to the 
	molar concentration, we have to convert them into the ones related to
        mass concentration first.  dc = (kt/T)dT  */

        ktmixmass = ktmix*dcdx;
        ktmpmass = ktmp*dcdx;

	psi = - betamass*ktmixmass/(alp*tk);  
	psip = - betamass*ktmpmass/(alp*tk);

   /*  factor of Rayleigh number alpha*g/(kapa*nu)   */

        Rf = alp*(G/100.)/(kapap*nu);
	s_t = 1.0/(Rf*pow(thick, 3));
	deltaT = 1707.76*s_t;

	s_c = kapap*nu/(betamass*(G/100.)*pow(thick, 3));

   /*  Dufour number  */  

	Q = tk*alp*alp/(betamass*betamass);
	Q = Q*dmudcm*(x_h*MW_HE + x_sf*MW_SF6)/(cpmix*1000.0);    

/* Output */

	*pro = 1.e-3*ro;		/* density in g/cm^3 */
	*palp = alp;			/* expansion coeff in 1/K */
	*pcomp = 0.;			/* compressibility, not calculated here */
	*plambmix = 1.e5*lambmix;	/* conductivity, in erg/s cm K */
	*petamix = 10.*etamix;
	*pcpmix = 1.e7*cpmix;
	*pkapap = *plambmix/((*pro)*(*pcpmix));	/* thermal diffusivity in cm^2/s */

	return;
}

/* ************************************************************************ 
   This subroutine calculates correlation functions (collision integrals) 	
	omega_1 and omega_2.

    Temperature T > 208 K.
  ************************************************************************ */
omega(tstar, rostar, vstar, ome_1, ome_2)
double tstar, rostar, vstar, *ome_1, *ome_2;
  {
	double om_1, om_2;   		/*  correlation functions defined 
						in kinetic theory  */
	double alpha, a_2, a_3, a_4, tlo, al_10, ar, arm2;
	double b_2, b_4, b_6;
/*	double exp(), pow(), log();   */

	tlo = log(tstar);
	alpha = log(vstar) - log(tstar);
	al_10 = log(vstar) - log(10);

/*  for He and SF6-He tstar > 10 for T > 208 K  */

	if(tstar <= 10)
	  {       
	    om_1 = exp(0.295402 - 0.510069*tlo + 0.189395*tlo*tlo - 0.045427*pow(tlo, 3) + 0.0037928*pow(tlo, 4));
	    om_2 = exp(0.46641 - 0.56991*tlo + 0.19591*tlo*tlo - 0.03879*pow(tlo,3) + 0.00259*pow(tlo,4));
	  }
	else
	  {
	    ar = al_10*rostar;
	    arm2 = 1/(ar*ar);
	    a_2 = -33.0838 + arm2*(20.0862 + (72.1059/al_10) + pow(8.27648/al_10, 2));
	    a_3 = 101.571 - arm2*(56.4472 + (286.393/al_10) + pow(17.7610/al_10, 2));
	    a_4 = -87.7036 + arm2*(46.3130 + (277.146/al_10) + pow(19.0573/al_10, 2));

	    om_2 = (rostar*rostar)*(alpha*alpha)*(1.04 + a_2/(tlo*tlo) + a_3/pow(tlo, 3) + a_4/pow(tlo, 4));

	    b_2 = -267.0 + arm2*(201.57 + (174.672/al_10) + pow(7.36916/al_10, 2));
	    b_4 = 26700.0 - arm2*(19.2265 + (27.6938/al_10) + pow(3.29559/al_10, 2))*1000.0;
	    b_6 = -8.9e+5 + arm2*(6.31013 + (10.2266/al_10) + pow(2.33033/al_10, 2))*1.0e+5;

	    om_1 = (rostar*rostar)*(alpha*alpha)*(0.89 + b_2/(tstar*tstar) + b_4/pow(tstar, 4) + b_6/pow(tstar, 6));
	  }

	*ome_1 = om_1;
	*ome_2 = om_2;
	return;
  }	

/* *********************************************************************** 
   This subroutine calculates functions related to second virial coefficient.
 ************************************************************************** */
bstar(tstar, rostar, vstar, bz)
double tstar, rostar, vstar, *bz;
  {
	double c_0=1.21852, c_1=-1.91794, c_2=0.92850, c_3=-0.16763, c_4=0.00673;
	double alpha, beta, tlo, al_10, be_10, ag_10;
	double d_2, d_4, d_6, dm, rm;
/*	double exp(), pow(), sqrt(), log();   */

	tlo = log(tstar);
	alpha = log(vstar) - log(tstar);
	al_10 = log(vstar) - log(10);
	ag_10 = al_10 + GAMA;
	beta = pow(rostar, 3)*(pow(alpha+GAMA, 3) + (PI*PI/2.0)*(alpha + GAMA) + 2.40411);
	be_10 = pow(rostar, 3)*(pow(al_10+GAMA, 3) + (PI*PI/2.0)*(al_10 + GAMA) + 2.40411); 

/*	fprintf(stderr, "%lf \n", tstar);  */
	if(tstar <= 10)
	  {
	   *bz = -sqrt(tstar)*exp(1.0/tstar)*(c_0 +c_1*tlo + c_2*tlo*tlo + c_3*pow(tlo, 3) + c_4*pow(tlo, 4));
/*	   fprintf(stderr, "%lf \n", *bz);	*/
	  }
	else
	 {
	   rm = pow(rostar, 3)/(be_10*be_10);
 	   dm = pow(rostar*rostar/be_10, 3)*pow(3*ag_10*ag_10 + PI*PI/2.0, 2);
	   d_2 = -15.9057 + 9.85958/be_10 + rm*(25.6607*ag_10*ag_10 - 9.73766*ag_10 + 42.2102) + 3.24589*dm;
	   d_4 = 84.3304 -61.9124/be_10 + rm*(-227.258*ag_10*ag_10 + 103.256*ag_10 - 373.824) - 34.4187*dm;
	   d_6 = -149.037 + 119.937/be_10 + rm*(483.571*ag_10*ag_10 - 273.727*ag_10 + 795.442) + 91.2423*dm;

	   *bz = beta*(1 + d_2/(tlo*tlo) + d_4/pow(tlo, 4) + d_6/pow(tlo, 6));
/* 	   fprintf(stderr, "L = %lf \n", *bz);   */
	 }

	return;
  }

/* ********************************************************************** 
   Calculation of the second virial coefficient of the mixture. 
 ************************************************************************ */

virial(x_sf, ttk, b_mix, c_mix, d_mix)
double x_sf, ttk, *b_mix, *c_mix, *d_mix;
  {
    	double x_h, b_h, b_sf, b_hsf, b_m, tk;
	double c_sf, c_h, d_sf, c_m, d_m;
	double bstar_h, bstar_sf, bstar_hsf;
    	double tstar_h, tstar_sf, tstar_hsf;
/*	double pow();   */

    	x_h = 1.0 - x_sf;
	tk = ttk;
	tstar_h = ttk/EK_HE;
	tstar_sf = ttk/EK_SF6;
	tstar_hsf = ttk/EK_HSF6;

	bstar(tstar_h, ROSTAR_HE, VSTAR_HE, &bstar_h);
	bstar(tstar_sf, ROSTAR_SF6, VSTAR_SF6, &bstar_sf);
	bstar(tstar_hsf, ROSTAR_HSF6, VSTAR_HSF6, &bstar_hsf);

	b_h = (2.0*PI/3.0)*NA*pow(SIGMA_HE, 3)*bstar_h*1000.0;
/*	fprintf(stderr, "b_h = %lf\n", b_h);   */
	b_sf = (2.0*PI/3.0)*NA*pow(SIGMA_SF6, 3)*bstar_sf*1000.0;
/*	fprintf(stderr, "b_sf = %lf\n", b_sf);   */
	b_hsf = (2.0*PI/3.0)*NA*pow(SIGMA_HSF6, 3)*bstar_hsf*1000.0;

/*  !!! BUT the calculated b_hsf is significantly larger than measurements,
    I correct it here for the range of temperature used in our experiments. */

        b_hsf = b_hsf - 0.018;
/*	fprintf(stderr, "b_hsf = %lf\n", b_hsf);   */

	b_m = x_h*x_h*b_h + x_sf*x_sf*b_sf + 2.0*x_h*x_sf*b_hsf;

/* the following C and D are based on pure SF6 and pure He, see programs
helium.c and sf6.c */

	c_h = -11.6902030 + (4.04273550e+4)/tk - (1.93679827e+6)/(tk*tk);
        c_h = 1.0e-6*c_h;

	c_sf = 	MW_SF6*MW_SF6*MW_SF6*(-1.64315e-4 + 0.100787/tk - (3.03945e+3)/pow(tk,3))/R;
	d_sf = 	MW_SF6*MW_SF6*MW_SF6*MW_SF6*(6.674565e-7 - 3.714493e-4/tk +15.73237/pow(tk,3))/R;

	c_m = c_h*x_h + c_sf*pow(x_sf, 2.6);    /* follow hc.c, I here chose
						2.6 instead of 2.389 there to 	
						include weaker effects of c_112 
						and c_211, since SF6 is a 
						sphereical molecule  */ 
	d_m = d_sf*pow(x_sf, 4.0);

	*b_mix = b_m;
	*c_mix = c_m;
	*d_mix = d_m;

	return;
  }

/* ********************************************************************* 
   dynamic viscosity, thermal conductivity and heat capacity.
 *********************************************************************** */
vislamb(x_sf, ttk, ppp, bmix, eta_mix, lamb_mix, cpmix)
double x_sf, ttk, ppp, bmix, *eta_mix, *lamb_mix, *cpmix;
  {
    	double x_h, H_h, H_sf, H_hsf;
	double astar_hsf, mu_hsf, fm_h, fm_sf, estar_h, estar_sf;
    	double tstar_h, tstar_sf, tstar_hsf, mu_h, mu_sf;
	double dt=0.001, dtstar, ddtstar; 
	double om_1, om_2, dom_2, dom_1_hsf, dom_2_hsf, ddom_1_hsf, ddom_2_hsf;
	double astar_h, astar_sf, cstar_hsf, bstar_hsf, bstar_hsf2;
	double la_mon_sf, la_mon_hsf, la0_sf, la0_h, tt, theta;
	double mrd, tom, mtop, mbottom, B_sf, ppsf, pp, cp_sf_t, cp_sf_tp, cp_h;
	double la_mon_mix, Z, X, Y, la_int_mix, amix, vo;
	double om_1_hsf, om_2_hsf;
	double umid, u_h, u_sf, muyz, uy, uz, d_hsf, d_sf, c_mix;
/*	double pow(), log(), sqrt();   */

	x_h = 1.0 - x_sf;
	tstar_h = ttk/EK_HE;
	tstar_sf = ttk/EK_SF6;
	tstar_hsf = ttk/EK_HSF6;

/*  Viscosity at zero-density is based on the kinetic theory by J. Bzowski et al, J. Phys. Chem. Ref.  Data, V.19, 1179 (1990).  */

 /* viscosity of helium */
	omega(tstar_h, ROSTAR_HE, VSTAR_HE, &om_1, &om_2);
	astar_h = om_2/om_1;
	dtstar = tstar_h + dt;
	omega(dtstar, ROSTAR_HE, VSTAR_HE, &om_1, &dom_2);

	estar_h = 1.0 + (tstar_h/4.0)*(log(dom_2) - log(om_2))/dt;
	fm_h = 1.0 + (3.0/196.0)*pow(8.0*estar_h - 7.0, 2);
	mu_h = (5.0/16.0)*fm_h*sqrt(MW_HE*K*ttk/(PI*NA*1000.0))/(SIGMA_HE*SIGMA_HE*om_2);  

 /* viscosity of sf6 */
	omega(tstar_sf, ROSTAR_SF6, VSTAR_SF6, &om_1, &om_2);
	astar_sf = om_2/om_1;
	dtstar = tstar_sf + dt;
	omega(dtstar, ROSTAR_SF6, VSTAR_SF6, &om_1, &dom_2);
  
	estar_sf = 1.0 + (tstar_sf/4.0)*(log(dom_2) - log(om_2))/dt;
	fm_sf = 1.0 + (3.0/196.0)*pow(8.0*estar_sf - 7.0, 2);
	mu_sf = (5.0/16.0)*fm_sf*sqrt(MW_SF6*K*ttk/(PI*NA*1000.0))/(SIGMA_SF6*SIGMA_SF6*om_2);   

 /* mixture viscosity */

	mrd = MW_HE*MW_SF6/(MW_HE + MW_SF6);
	tom = MW_HE + MW_SF6;
	omega(tstar_hsf, ROSTAR_HSF6, VSTAR_HSF6, &om_1_hsf, &om_2_hsf);
	astar_hsf = om_2_hsf/om_1_hsf;

	dtstar = tstar_hsf + DDT;
	omega(dtstar, ROSTAR_HSF6, VSTAR_HSF6, &dom_1_hsf, &dom_2_hsf);
	ddtstar = tstar_hsf - DDT;
	omega(ddtstar, ROSTAR_HSF6, VSTAR_HSF6, &ddom_1_hsf, &ddom_2_hsf);

	cstar_hsf = 1.0 + (tstar_hsf/3.0)*(log(dom_1_hsf) - log(om_1_hsf))/DDT;
	
	bstar_hsf = 1.0 + 3.0*cstar_hsf - 3.0*cstar_hsf*cstar_hsf;
	bstar_hsf2 = log(dom_1_hsf) + log(ddom_1_hsf) - 2.0*log(om_1_hsf);
	bstar_hsf = bstar_hsf - tstar_hsf*tstar_hsf*bstar_hsf2/(3.0*DDT*DDT);

	mu_hsf = (5.0/16.0)*sqrt(2.0*mrd*K*ttk/(PI*NA*1000.0))/(SIGMA_HSF6*SIGMA_HSF6*om_2_hsf);
	H_h = x_h*x_h/mu_h + (2.0*x_h*x_sf/mu_hsf)*(MW_HE*MW_SF6/(tom*tom))*(5.0/(3.0*astar_hsf) + MW_SF6/MW_HE);
	H_sf = x_sf*x_sf/mu_sf + (2.0*x_h*x_sf/mu_hsf)*(MW_HE*MW_SF6/(tom*tom))*(5.0/(3.0*astar_hsf) + MW_HE/MW_SF6);
	H_hsf = - (2.0*x_h*x_sf/mu_hsf)*(MW_HE*MW_SF6/(tom*tom))*(5.0/(3.0*astar_hsf) - 1.0);

	mtop = 2.0*x_h*x_sf*H_hsf - x_h*x_h*H_sf - x_sf*x_sf*H_h;
	mbottom = H_h*H_sf - H_hsf*H_hsf;
	if(x_h == 0.0) *eta_mix = mu_sf;
        else if(x_h == 1.0) *eta_mix = mu_h;
	else *eta_mix = - mtop/mbottom;

/*   Thermal conductivity of mixture:  based on a model proposed for polyatomic gases by J. Kestin and W. A. Wakeham, Int. J. Thermophys. V.4, 295 (1983). They employed measured physical properties whenever possible.  A. A. Clifford, J. Kestin and W. A. Wakeham, Ber. Bunsenges. Phys. Chem. V.85, 385 (1981) also showed the method to be better than direct calculations from the kinetic theory.  However, the astar and bastar are still from the theory.    */

	tt = ttk - ZERO_C;

 /*  thermal conductivity of helium at zero-density. The formulas are copied from helium.c. */

	theta = tt/320.0;
  	la0_h = 0.159*pow(theta+0.8536, 0.69);
	la0_h = la0_h + 0.880*theta/(pow(theta, 5) + 128.0);

 /*  thermal conductivity of pure SF6 at zero-density. The formulas are copied from sf6.c.  */

	la0_sf = 0.01303*(1.0 + 0.00549*(tt - 27.5)); 

  /* 1--> calculate la_mon_mix, treating the polyatomic molecule as a monatomic one */
 
	la_mon_sf = 15.0*R*mu_sf/(4.0*MW_SF6);
	la_mon_hsf = 15.0*R*mu_hsf*(MW_HE + MW_SF6)/(8.0*MW_HE*MW_SF6);	

	umid = 4.0*astar_hsf/15.0 + 0.5*(MW_HE - MW_SF6)*(MW_HE - MW_SF6)/(MW_HE*MW_SF6);
	u_h = umid -(bstar_hsf/5.0 + 1.0/12.0)*MW_HE/MW_SF6;
        u_sf = umid -(bstar_hsf/5.0 + 1.0/12.0)*MW_SF6/MW_HE;

        muyz = (MW_SF6 + MW_HE)*(MW_SF6 + MW_HE)/(4.0*MW_SF6*MW_HE);

	uy = (4.0*astar_hsf/15.0)*muyz*la_mon_hsf*la_mon_hsf/(la_mon_sf*la0_h) - bstar_hsf/5.0 - 1.0/12.0;
        uy = uy - (12.0*bstar_hsf - 25.0)*(MW_HE - MW_SF6)*(MW_HE - MW_SF6)/(32.0*astar_hsf*MW_HE*MW_SF6);

        uz = (4.0*astar_hsf/15.0)*(muyz*(la_mon_hsf/la_mon_sf + la_mon_hsf/la0_h) - 1.0) - bstar_hsf/5.0 - 1.0/12.0;

        X = x_h*x_h/la0_h + 2.0*x_h*x_sf/la_mon_hsf + x_sf*x_sf/la_mon_sf;
        Y = x_h*x_h*u_h/la0_h + 2.0*x_h*x_sf*uy/la_mon_hsf + x_sf*x_sf*u_sf/la_mon_sf;
	Z =  x_h*x_h*u_h + 2.0*x_h*x_sf*uz + x_sf*x_sf*u_sf;

	la_mon_mix = (1.0 + Z)/(X + Y);

  /* 2 ---> find out the diffusion constants for internal energy in order to calculate the contribution from internal degree of freedom.  */

	d_hsf = 3.0*astar_hsf*(MW_SF6 + MW_HE)*mu_hsf/(5.0*MW_SF6*MW_HE);
	d_sf = 6.0*astar_sf*mu_sf/(5.0*MW_SF6);

	/* fprintf(stderr, " d_hsf = %.3e  d_sf = %.3e\n", d_hsf, d_sf);  */
	
	la_int_mix = (la0_sf - la_mon_sf)/(1.0 + x_h*d_sf/(x_sf*d_hsf));

  /* 3 ---> final thermal conductivity */

	*lamb_mix = la_mon_mix + la_int_mix;

/* Heat capacity */

/*  the volumn of 1 mole mixture  */

	amix = sqrt(1. + 4.*ppp*bmix/(R*ttk));
	vo = 0.5*R*ttk*(1.+amix)/ppp;     

 /*  partial pressure of SF6. B_sf is copied from sf6-p.c.  */

        B_sf = MW_SF6*MW_SF6*(0.054820 - 35.874/ttk - 1.1443e+6/pow(ttk,3))/R;
	ppsf = x_sf*R*ttk*(1.0 + x_sf*B_sf/vo)/vo;
	pp = ppsf/1.0e+5;	             /*  unit: bar  */

 /*  heat capacity.  unit:  kJ/kmole*K.  The formulas for the isobaric heat
       capacity of SF6 are copied from sf6-p.c.  */
	
	cp_h = 20.78585; 	   /* pure helium isobaric heat capacity */

	cp_sf_t = 91.07 + 0.2581*tt - 4.504e-4*tt*tt;
	cp_sf_tp = cp_sf_t*(1.0 + (0.0023 + 1.346e-5*tt)*(pp-1.0) + (63.56 - 0.7144*tt)*1.0e-5*(pp-1.0)*(pp-1.0));

  /*  heat capacity of the mixture  */

	c_mix = cp_h*x_h + cp_sf_tp*x_sf;

	*cpmix = c_mix;

	return;
}

/* ******************************************************************** 
   Calculation of the diffusion coefficients for pure Helium, SF6 and their
   mixture.  This subroutine also calculates the thermal diffusion ratio. 

   The self-diffusion constants are going to be used to calculate thermal 
   conductivity of mixture.
 ******************************************************************** */
diff(x_sf, ttk, ppp, D_h, D_sf, D_mix, rd_ktmix)
double x_sf, ttk, ppp, *D_h, *D_sf, *D_mix, *rd_ktmix;
  { 
	double x_h, delta_hsf, kt, dh, dsf;
	double astar_hsf, bstar_hsf, cstar_hsf, a_hsf, b_hsf, a_m, m_sf2h, m_m;
	double fd_h, fd_sf, astar_h, cstar_h, astar_sf, cstar_sf;
    	double tstar_h, tstar_sf, tstar_hsf, bstar_hsf2;
	double dtstar, ddtstar, dtstar_h, dtstar_sf;
	double s_h, s_sf, q_h, q_sf, sq_m_sf, sq_m_h, q_hsf1, q_hsf2, q_hsf; 
	double om_1_hsf, om_2_hsf, dom_1_hsf, dom_2_hsf;
	double ddom_1_hsf, ddom_2_hsf, om_1_h, om_2_h, om_1_sf, om_2_sf;
	double dom_1_h, dom_2_h, dom_1_sf, dom_2_sf;
	double dm, mplus, rkt;
/*	double pow(), log(), sqrt();   */

	m_sf2h = MW_HE/MW_SF6;
	m_m = MW_HE*MW_SF6/(MW_HE + MW_SF6);
	mplus = MW_HE + MW_SF6;

	tstar_hsf = ttk/EK_HSF6;
	omega(tstar_hsf, ROSTAR_HSF6, VSTAR_HSF6, &om_1_hsf, &om_2_hsf);
	dtstar = tstar_hsf + DDT;
	omega(dtstar, ROSTAR_HSF6, VSTAR_HSF6, &dom_1_hsf, &dom_2_hsf);
	ddtstar = tstar_hsf - DDT;
	omega(ddtstar, ROSTAR_HSF6, VSTAR_HSF6, &ddom_1_hsf, &ddom_2_hsf);
	/*  fprintf(stderr, "om_1_hsf= %f \n", om_1_hsf);  */

        tstar_h = ttk/EK_HE;
	omega(tstar_h, ROSTAR_HE, VSTAR_HE, &om_1_h, &om_2_h);
	dtstar_h = tstar_h + DDT;
	omega(dtstar_h, ROSTAR_HE, VSTAR_HE, &dom_1_h, &dom_2_h);

	tstar_sf = ttk/EK_SF6;
	omega(tstar_sf, ROSTAR_SF6, VSTAR_SF6, &om_1_sf, &om_2_sf);
	dtstar_sf = tstar_sf + DDT;
	omega(dtstar_sf, ROSTAR_SF6, VSTAR_SF6, &dom_1_sf, &dom_2_sf);

	a_m = (1.8*m_sf2h + 1.0)*(1.8*m_sf2h + 1.0)*8.0*om_2_h;
	a_hsf = sqrt(2.0)*om_1_hsf/a_m;
	b_hsf = 10.0*a_hsf*(1.0 + 1.8*m_sf2h + 3*m_sf2h*m_sf2h) - 1.0;
	
	cstar_hsf = 1.0 + (tstar_hsf/3.0)*(log(dom_1_hsf) - log(om_1_hsf))/DDT;
	/*  fprintf(stderr, "cstar_hsf= %f \n", cstar_hsf);  */
	
	bstar_hsf = 1.0 + 3.0*cstar_hsf - 3.0*cstar_hsf*cstar_hsf;
	bstar_hsf2 = log(dom_1_hsf) + log(ddom_1_hsf) - 2.0*log(om_1_hsf);
	bstar_hsf = bstar_hsf - tstar_hsf*tstar_hsf*bstar_hsf2/(3.0*DDT*DDT);

	delta_hsf = 1.3*(6.0*cstar_hsf - 5.0)*(6.0*cstar_hsf - 5.0);
	delta_hsf = delta_hsf*a_hsf*x_sf/(1 + b_hsf*x_sf);
/*	fprintf(stderr, "delta_hsf= %f \n", delta_hsf);   */

	dm = (1.0 + delta_hsf)/(SIGMA_HSF6*SIGMA_HSF6*om_1_hsf);
	dm = dm*(3.0/8.0)*K*ttk/ppp;
        dm = dm*sqrt(K*ttk/(2.0*PI*m_m/(NA*1000.0)));

	sq_m_sf = sqrt(m_m*2.0/MW_SF6)*om_2_sf*SIGMA_SF6*SIGMA_SF6;
	sq_m_sf = sq_m_sf/(SIGMA_HSF6*SIGMA_HSF6*om_1_hsf);

	sq_m_h = sqrt(m_m*2.0/MW_HE)*om_2_h*SIGMA_HE*SIGMA_HE;
	sq_m_h = sq_m_h/(SIGMA_HSF6*SIGMA_HSF6*om_1_hsf);

	astar_hsf = om_2_hsf/om_1_hsf; 

	s_sf = (MW_SF6/MW_HE)*sq_m_sf - 4.0*m_m*astar_hsf/mplus;
	s_sf = s_sf + 7.5*MW_HE*(MW_SF6 - MW_HE)/(mplus*mplus);

	s_h = (MW_HE/MW_SF6)*sq_m_h - 4.0*m_m*astar_hsf/mplus;
	s_h = s_h + 7.5*MW_SF6*(MW_HE - MW_SF6)/(mplus*mplus);

/*	fprintf(stderr, "s_h= %f \n", s_h);
	fprintf(stderr, "s_sf= %f \n", s_sf);  */

	q_sf = 2.0*sq_m_sf/(MW_HE*mplus);
	q_sf = q_sf*((2.5 - 1.2*bstar_hsf)*MW_SF6*MW_SF6 + 3.0*MW_HE*MW_HE + 1.6*MW_SF6*MW_HE*astar_hsf);

	q_h = 2.0*sq_m_h/(MW_SF6*mplus);
	q_h = q_h*((2.5 - 1.2*bstar_hsf)*MW_HE*MW_HE + 3.0*MW_SF6*MW_SF6 + 1.6*MW_SF6*MW_HE*astar_hsf);

	q_hsf1 = 15.0*(2.5 - 1.2*bstar_hsf)*(MW_HE - MW_SF6)*(MW_HE - MW_SF6);
	q_hsf1 = q_hsf1/(mplus*mplus);
	q_hsf1 = q_hsf1 + 4.0*m_m*astar_hsf*(11.0 - 2.4*bstar_hsf)/mplus;
	q_hsf2 = 1.6*sqrt(mplus/m_m)*SIGMA_HE*SIGMA_HE*om_2_h;
	q_hsf2 = q_hsf2*SIGMA_SF6*SIGMA_SF6*om_2_sf/pow(SIGMA_HSF6, 4);
	q_hsf2 = q_hsf2/(om_1_hsf*om_1_hsf);
	q_hsf = q_hsf1 + q_hsf2;

	x_h = 1.0 - x_sf;

	kt = (6.0*cstar_hsf - 5.0)*(x_sf*s_sf - x_h*s_h);
	rkt = kt/(x_sf*x_sf*q_sf + x_h*x_h*q_h + x_sf*x_h*q_hsf);
/*	fprintf(stderr, "s_sf=%.3e  s_h=%.3e\n", s_sf, s_h);    */

	*rd_ktmix = rkt;	    /*  reduced thermal diffusion factor */
	*D_mix = dm;

/*  now the self-diffusion constants for calculating thermal conductivity
    of mixture using Cheung et al's method.  Not useful anymore, but I still
    leave them here.    */

	astar_h = om_2_h/om_1_h;
	cstar_h = 1.0 + (tstar_h/3.0)*(log(dom_1_h) - log(om_1_h))/DDT; 

	astar_sf = om_2_sf/om_1_sf;
	cstar_sf = 1.0 + (tstar_sf/3.0)*(log(dom_1_sf) - log(om_1_sf))/DDT;

	fd_h = (6.0*cstar_h - 5.0)*(6.0*cstar_h - 5.0)/(2.0*astar_h + 5.0);
	fd_h = 1.0 + fd_h/8.0;

	fd_sf = (6.0*cstar_sf - 5.0)*(6.0*cstar_sf - 5.0)/(2.0*astar_sf + 5.0);
	fd_sf = 1.0 + fd_sf/8.0;

	dh = (3.0/8.0)*K*ttk*fd_h/ppp;
        dh = dh*sqrt(K*ttk/(PI*MW_HE/(NA*1000.0)))/(SIGMA_HE*SIGMA_HE*om_1_h);

	dsf = (3.0/8.0)*K*ttk*fd_sf/ppp;
        dsf = dsf*sqrt(K*ttk/(PI*MW_SF6/(NA*1000.0)))/(SIGMA_SF6*SIGMA_SF6*om_1_sf);

	/* fprintf(stderr, " dh = %.3e  dsf = %.3e\n", dh, dsf);   */

	*D_h = dh;
	*D_sf = dsf;
	
	return;
}

/* ***********************************************************************
   Calculation of the mixture density, thermal expansion, solutal expansion
   and baric expansion coefficients, and refractive index.
 ********************************************************************** */
rhoalp(x_sf, tk, ppa, rho, alp, beta, ksi, n, n_c_hat, n_t_hat)
double x_sf, tk, ppa, *rho, *alp, *beta, *ksi, *n, *n_c_hat, *n_t_hat;
  {
	double ro, al, x_h, dtk, dro, bet, dx_sf;
	double dppa, ks;
	double cmass_sf, rho1, rho2, dcmass_sf, drho1, drho2; 
	double ind, dind, middle, middle2, mixmiddle, nn, nc, nt;

	double super_vir();

	x_h = 1.0 - x_sf;
	
	ro = super_vir(x_sf, tk, ppa);				  /* kg/m^3 */

	cmass_sf = x_sf*MW_SF6/(x_sf*MW_SF6 + x_h*MW_HE); 
						       /* mass concentration */
	rho1 = cmass_sf*ro*1000.0/MW_SF6;  /* molar density of SF6, mole/m^3 */ 
	rho2 = (1.0-cmass_sf)*ro*1000.0/MW_HE;	     /*  molar density of He */  
        middle = AR_SF6*rho1 + BR_SF6*rho1*rho1;
        middle2 = AR_HE*rho2 + BR_HE*rho2*rho2;
	mixmiddle = middle + middle2;
        ind = sqrt(3.0/(1.0-mixmiddle) - 2.0);
	nn = ind;				      /*  refractive index */

/* derivatives with respect to temperature */

	dtk = tk + DDT;
	dro = super_vir(x_sf, dtk, ppa);

	al = (ro - dro)/(DDT*ro);                      /* thermal expansion  */

	drho1 = cmass_sf*dro*1000.0/MW_SF6;         /* molar density of SF6  */ 
	drho2 = (1.0-cmass_sf)*dro*1000.0/MW_HE;     /*  molar density of He */  
        middle = AR_SF6*drho1 + BR_SF6*drho1*drho1;
        middle2 = AR_HE*drho2 + BR_HE*drho2*drho2;
	mixmiddle = middle + middle2;
        dind = sqrt(3.0/(1.0-mixmiddle) - 2.0);
	nt = (dind - ind)/DDT;

/* derivatives with respect to concentration */

	dx_sf = x_sf + DDC;
	dro = super_vir(dx_sf, tk, ppa);

	bet = (ro - dro)/(DDC*ro);    /* solutal expansion  */
	
	dcmass_sf = dx_sf*MW_SF6/(dx_sf*MW_SF6 + (1.-dx_sf)*MW_HE); 
						       /* mass concentration */
	drho1 = dcmass_sf*dro*1000.0/MW_SF6;        /* molar density of SF6  */
	drho2 = (1.0-dcmass_sf)*dro*1000.0/MW_HE;    /*  molar density of He */   
        middle = AR_SF6*drho1 + BR_SF6*drho1*drho1;
        middle2 = AR_HE*drho2 + BR_HE*drho2*drho2;
        mixmiddle = middle + middle2;
        dind = sqrt(3.0/(1.0-mixmiddle) - 2.0);
	nc = (dind - ind)/(dcmass_sf - cmass_sf);

/* derivatives with respect to pressure */

	dppa = ppa + DDP;
	dro = super_vir(x_sf, tk, dppa);
	ks = (dro - ro)/(DDP*ro);    /* pressure contraction  */

	*rho = ro;
	*alp = al;
	*beta = bet;
	*ksi = ks;
	*n = nn;
	*n_c_hat = nc;
	*n_t_hat = nt;

	return;
  }

/* **************************************************************************
  The following subroutine calculates the pressure dependence of diffusion
  coefficient Dmix and thermal diffusion factor ktmix.  These calculations
  are based on the measurements of T. N. Bell et al. (Chem. Phys. Lett. 
  V.45, 445 (1977)) and R.D. Trengove et al. (Physica A V.108, 502 (1981)). The 
  experiments only went up to 9 atm for Dmix and 5 atm for ktmix.  However, I   
  generalize the "second virial coefficients" from these experiments to our 
  higher pressures.  
 ************************************************************************* */
dktpress_he_sf6(ppa, x_sf, ro, ttk, Dmix, rkt, Dp, kmix, ktp, rckt, rktp)
double ppa, x_sf, ro, ttk, Dmix, rkt, *Dp, *kmix, *ktp, *rckt, *rktp;
{
	double dmp, patm, x_h;
	double bkt, kt, nd, rct, rpt;
	double bh=4.67e-4, bsf=-121.1e-4, be=72.94e-4; /* atm^-1 */
	double sh=0.193, ssf=0.650;   /* nm */
	double ff1, ff2;

  	patm = ppa/101330.0;
	x_h = 1.0 - x_sf;

	nd = ro/(MW_SF6*x_sf + MW_HE*(1.0 - x_sf));
	nd = nd*NA*1000.0;			          /* unit: 1/m^3   */

  /*  Dmix ~ p, T.N. Bell, I. R. Shankland and P. J. Dunlop, Chem. Phys. Lett.
      V.45, 445 (1977). */

	ff1 = 1.0 - 4.0*x_sf*x_h*be*patm;
	ff2 = (x_sf*bsf*(4.0*sh + ssf) + x_h*bh*(4.0*ssf + sh))/((sh + ssf)*101330.0);
	ff2 = 1.0 + 0.25*nd*K*ttk*ff2;
/*	printf("ff1=%f\n", ff1);  */
/*	printf("ff2=%f\n", ff2);  */
	dmp = Dmix*ff1/ff2;
	*Dp = dmp;

  /*  ktmix ~ pressure (R.D. Trengove et al. Physica A V.108, 502(1981)).  Here
      we also correct the error of theoretical calculations which are about 10%
      smaller than available data.  */

	rct = rkt*(1.1106 + 0.06796*x_sf);    /* fit and correct the difference 
							between experiments and
							calculations   */
	*rckt = rct;	      /* corrected reduced thermal diffusion factor */

	kt = rct*x_sf*(1.0 - x_sf);	/* not reduced */
	*kmix = kt;

	bkt =  1.9552e-27*x_sf  +  3.9344e-27*x_sf*x_sf;   
		    /*   fit the expermients of Trengove et al. at 300 K   */	
	kt = kt*(1.0 + bkt*nd);
	rpt = rct*(1.0 + bkt*nd);        /*  pressure-dependent of k_t  */
	
        *ktp = kt;
	*rktp = rpt;

        return;
}
	
/* ****************************************************************************
   This subroutine calculates the derivative of the chemical potential of an 
   binary gas mixture to the molar concentration.  Since I found no 
   experimental data, we use the formula from the ideal gas mixture (See
   Landau & Lifishitz, Fluid Mechanics, (Pergamon, 1959), p221 and p226).
   In Dufour number, the derivative is respect to mass concentration (see Hort 
   et al. PRA V.45, 3737 (1992)). 
 *************************************************************************** */
chemical_he_sf6(tk, x_sf, dmudcmass)
double tk, x_sf, *dmudcmass;
  {
	double cmass_sf, dmdc;    /* dmdc = d(mu)/dc, mu is the chemical 
						potential  */
	
	cmass_sf = x_sf*MW_SF6/(x_sf*MW_SF6 + (1.0 - x_sf)*MW_HE);
	dmdc = R*tk/((1.0 - cmass_sf)*MW_SF6 + cmass_sf*MW_HE);
	dmdc = dmdc*(1.0/cmass_sf + 1.0/(1.0 - cmass_sf));
					
/*     fprintf(stderr, "cmass=%f dmudcmass = %.4e J/kg \n", cmass_sf, dmdc);   */					
	*dmudcmass = dmdc;

  	return;
  }	

/* **************************************************************************
   subroutines solving the equation of state.
*************************************************************************** */

double super_vir(x_sf,temp,press)
double x_sf, temp, press;
  {
	double rho, mw;
	double rho0;
	double xacc,x1,x2;	
	double rtnewt();

	xacc = 1.e-5;
        mw = x_sf*MW_SF6 + (1.0-x_sf)*MW_HE;
	rho0 = press*mw/(R*temp);
	x1 = 0.1*rho0;
	x2 = 1.9*rho0;
	rho = rtnewt(x1,x2,xacc,x_sf,temp,press);
	return rho;
  }

double rtnewt(x1,x2,xacc,x_sf,temp,press)
double x1,x2,xacc,x_sf,temp,press;
{
	int j;
	double df,dx,f,rtn;

	rtn=0.5*(x1+x2);
	for (j=1;j<=JMAX;j++) {
		funcd(rtn,x_sf,temp,press,&f,&df);
		dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0){
			printf("Jumped out of brackets in RTNEWT\n");
			exit(0);
		}
		if (fabs(dx) < xacc) return rtn;
	}
	printf("Maximum number of iterations exceeded in RTNEWT\n");
}

funcd(rho,x_sf,temp,press,func,dfunc)
double rho, x_sf, temp, press, *func, *dfunc;
  {
	double b_mix, c_mix, mw, bb, cc, d_mix, dd;
	double fun, dfun;

	virial(x_sf, temp, &b_mix, &c_mix, &d_mix);
	mw = x_sf*MW_SF6 + (1.0-x_sf)*MW_HE;
	bb = R*temp*b_mix/(mw*mw);
        cc = R*temp*c_mix/pow(mw, 3.0);
	dd = R*temp*d_mix/pow(mw, 4.0);

	fun = R*temp*rho/mw + bb*rho*rho + cc*pow(rho, 3.0); 
        fun = fun + dd*pow(rho, 4.0) - press;
	dfun = R*temp/mw + 2.*bb*rho + 3.*cc*rho*rho;
        dfun = dfun + 4.*dd*pow(rho, 3.0);

	*func = fun;
	*dfunc = dfun;
	return;
}
