/*
 * main.c
 * Copyright (C) 2016  
 * 
 * bianchi_phi is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * bianchi_phi is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */  

#include <math.h>
#include <float.h>
#include <stdlib.h> 
#include <stdio.h>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <gsl/gsl_spline.h>

int cvode_check_flag (void *flagvalue, const char *funcname, int opt);
double sqrt1px_m1 (const double x);

#define LOG_DBL_MAX 7.0978271289338397e+02

#define CVODE_CHECK(chk,name,val,ret) \
do { \
  if (!cvode_check_flag (chk, name, val)) \
    return ret; \
} while (FALSE)

/********************************************************* 
 *********************** STRUCTURES **********************
*********************************************************/

struct _BianchiConsts
{
	double gamma;
	double beta;
	double g0;
	double bg;
	double V0;
	double bV;
	double p;
	double q;
	double sigma2_0;
};

typedef struct _BianchiConsts BianchiConsts;

struct _SyState 
{
	double a;
	double alpha;
	double phi;
	double varphi;
	double H;
	double bounce;
	int j_bounce;
	double alpha_i;
	double alpha_f;
	double t_0;
	int    i;
	int cont;
	double phi_inf;
	double phi_sup;
	double varphi_inf;
	double varphi_sup;
	//~ int    a_size;
};

typedef struct _SyState SyState;

/*******************************************************
*************** PRÉ-DEFINED FUCTIONS *******************
********************************************************/

static double Bianchi_P (const double t, const double a, const double phi, const double varphi, const double H, BianchiConsts *c1);
static double Bianchi_D (const double t, const double a, const double phi, const double varphi, const double H, BianchiConsts *c1);
//~ static double Bianchi_beta_p (const double eps, const double sign_beta_p, BianchiConsts *c1);


/*******************************************************
********************** FUCTIONS ************************
********************************************************/

static double
Bianchi_sigma2 (const double a, BianchiConsts *c1)
{
	const double a2 = a * a;
	const double a3 = a2 * a;
	const double a6 = a3 * a3;

	return c1->sigma2_0 / a6;
}

//~ static double 
//~ Bianchi_R(const double a, const double H, BianchiConsts *c1)
//~ {
	//~ const double H2			= H * H;
	//~ const double sigma_2	= Bianchi_sigma2(a, c1);
	//~ const double R			= H2 - sigma_2 / (6.0);
//~ 
	//~ return R;	
//~ }

//~ static double
//~ Bianchi_B( const double a, const double H, BianchiConsts *c1)
//~ {
	//~ const double R	= Bianchi_R(a, H, c1);
	//~ const double V0	= - c1-> V0;
	//~ const double B	=  (- 1.5 * R) / V0;
//~ 
	//~ return B;
//~ }

static double
Bianchi_Vg (const double phi, const double Vg0, const double b, const double pq)
{
	const double exp_arg1 = sqrt (2.0 / pq) * phi;
	const double exp_arg2 = b * exp_arg1;
	const double exp1     = exp ( - exp_arg1 );
	const double exp2     = exp ( + exp_arg2 );

	return 2.0 * Vg0 / (exp1 + exp2);
}

static double
Bianchi_dVg_dphi (const double phi, const double Vg0, const double b, const double pq)
{
	const double sqrt2p      = sqrt (2.0 / pq);
	const double exp_arg1    = sqrt2p * phi;
	const double exp_arg2    = b * exp_arg1;
	const double exp1        = exp ( - exp_arg1 );
	const double exp2        = exp ( + exp_arg2 );
	const double denom       = (exp1 + exp2);

	const double ddenom_dphi = - sqrt2p * exp1 + b * sqrt2p * exp2;

	if (-exp_arg1 > LOG_DBL_MAX)
	return 0.0;
	if (exp_arg2 > LOG_DBL_MAX)
	return 0.0;

	return - 2.0 * Vg0 * ddenom_dphi / (denom * denom);
}

static void
Bianchi_H (const double t, const double a, const double phi, const double varphi, BianchiConsts *c1, double *Hp, double *Hm)
{
	const double varphi2    = varphi * varphi;
	const double varphi3    = varphi2 * varphi;
	const double varphi4    = varphi2 * varphi2;
	const double sigma2     = Bianchi_sigma2 (a, c1);
	const double g          = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);
	const double V          = Bianchi_Vg (phi, c1->V0, c1->bV, c1->q);
	const double J          = 0.5 * (1.0 - g) * varphi2 + (3.0 / 4.0) * varphi4 * c1->beta + V;
	const double B          = - c1->gamma * varphi3;
	const double C          = - J * (1.0/ 3.0) - sigma2 * (1.0 / 6.0);
	const double Delta      = B * B - 4.0 * C;
	const double Delta_B2m1 = - 4.0 * C / (B * B);

	if (B > 0.0)
	{
		*Hp = B * sqrt1px_m1 (Delta_B2m1) * 0.5;
		*Hm = (- B - sqrt (Delta)) * 0.5;
	}
	else
	{
		*Hp = (- B + sqrt (Delta)) * 0.5;
		*Hm = B * sqrt1px_m1 (Delta_B2m1) * 0.5;
	}
}

static double
Bianchi_rho (const double t, const double a, const double phi, const double varphi, const double H, BianchiConsts *c1)
{
	const double varphi2 = varphi * varphi;
	const double varphi3 = varphi2 * varphi;
	const double varphi4 = varphi2 * varphi2;
	const double g       = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);
	const double V       = Bianchi_Vg (phi, c1->V0, c1->bV, c1->q);

	const double rho     = 0.5 * (1.0 - g) * varphi2 + (3.0 / 4.0) * c1->beta * varphi4 + 3.0 * c1->gamma * H * varphi3 + V;

	return rho;
}

static double
Bianchi_p (const double t, const double a, const double phi, const double varphi, const double H, BianchiConsts *c1)
{
	const double varphi2 = varphi * varphi;
	const double varphi4 = varphi2 * varphi2;
	const double g       = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);
	const double V       = Bianchi_Vg (phi, c1->V0, c1->bV, c1->q);
	const double dV      = Bianchi_dVg_dphi (phi, c1->V0, c1->bV, c1->q);
	const double P       = Bianchi_P (t, a, phi, varphi, H, c1);
	const double D       = Bianchi_D (t, a, phi, varphi, H, c1);
	const double ddphi   = - (D * varphi + dV) / P;

	const double p       = 0.5 * (1.0 - g) * varphi2 + 0.25 * c1->beta * varphi4 - c1->gamma * varphi2 * ddphi - V;

	return p;
}

static double
Bianchi_rho_p_p (const double t, const double a, const double phi, const double varphi, const double H, BianchiConsts *c1)
{
	const double varphi2 = varphi * varphi;
	const double varphi3 = varphi2 * varphi;
	const double varphi4 = varphi2 * varphi2;
	const double g       = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);
	const double dV      = Bianchi_dVg_dphi (phi, c1->V0, c1->bV, c1->q);
	const double P       = Bianchi_P (t, a, phi, varphi, H, c1);
	const double D       = Bianchi_D (t, a, phi, varphi, H, c1);
	const double ddphi   = - (D * varphi + dV) / P;

	const double rho_p_p = (1.0 - g) * varphi2 + c1->beta * varphi4 + 3.0 * c1->gamma * H * varphi3 - c1->gamma * varphi2 * ddphi;

	return rho_p_p;
}

static double
Bianchi_P (const double t, const double a, const double phi, const double varphi, const double H, BianchiConsts *c1)
{
	const double varphi2 = varphi * varphi;
	const double varphi4 = varphi2 * varphi2;
	const double gamma2  = c1->gamma * c1->gamma;
	const double g       = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);

	const double P       = (1.0 - g) + 6.0 * c1->gamma * H * varphi + 3.0 * c1->beta * varphi2 + 1.5 * gamma2 * varphi4;

	return P;
}

static double
Bianchi_D (const double t, const double a, const double phi, const double varphi, const double H, BianchiConsts *c1)
{
	const double varphi2 = varphi * varphi;
	const double varphi3 = varphi2 * varphi;
	const double varphi4 = varphi2 * varphi2;
	const double varphi5 = varphi3 * varphi2;
	const double gamma2  = c1->gamma * c1->gamma;
	const double sigma2 = Bianchi_sigma2 (a, c1);
	const double g      = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);
	const double dg     = Bianchi_dVg_dphi (phi, c1->g0, c1->bg, c1->p);
	const double H2     = H * H;

	const double D      = 3.0 * (1.0 - g) * H + (9.0 * c1->gamma * H2 - 0.5 * dg) * varphi
	+ 3.0 * c1->beta * H * varphi2 - 1.5 * (1.0 - g) * c1->gamma * varphi3 
	- 4.5 * gamma2 * H * varphi4 - 1.5 * c1->beta * c1->gamma * varphi5 
	- 1.5 * c1->gamma * sigma2 * varphi;

	return D;
}

static double
Bianchi_dH (const double t, const double a, const double phi, const double varphi, const double H, BianchiConsts *c1)
{
	const double sigma2  = Bianchi_sigma2 (a, c1);
	const double rho_p_p = Bianchi_rho_p_p (t, a, phi, varphi, H, c1);
	const double dH      = -0.5 * (rho_p_p + sigma2);

	return dH;
}

/*******************************************************
****************** RHS - cosmic time********************
********************************************************/

static int 
Bianchi_f_t(realtype t, N_Vector y_t, N_Vector ydot_t, void *f_data)
{
	BianchiConsts *c1    = (BianchiConsts *) f_data;
	const double a       = NV_Ith_S (y_t, 0);
	const double phi     = NV_Ith_S (y_t, 1);
	const double varphi  = NV_Ith_S (y_t, 2);
	const double H       = NV_Ith_S (y_t, 3);

	const double dV      = Bianchi_dVg_dphi (phi, c1->V0, c1->bV, c1->q);
	const double P       = Bianchi_P (t, a, phi, varphi, H, c1);
	const double D       = Bianchi_D (t, a, phi, varphi, H, c1);
	const double ddphi   = - (D * varphi + dV) / P;
	const double dH      = Bianchi_dH (t, a, phi, varphi, H, c1);

	NV_Ith_S (ydot_t, 0)   = H * a;
	NV_Ith_S (ydot_t, 1)   = varphi;
	NV_Ith_S (ydot_t, 2)   = ddphi;
	NV_Ith_S (ydot_t, 3)   = dH;

	return 0;
}

/*******************************************************
******************** RHS - alpha ***********************
********************************************************/

static int 
Bianchi_f_alpha(realtype alpha, N_Vector y_alpha, N_Vector ydot_alpha, void *f_data)
{
	BianchiConsts *c1    = (BianchiConsts *) f_data;
	const double a       = NV_Ith_S (y_alpha, 0);
	const double phi     = NV_Ith_S (y_alpha, 1);
	const double varphi  = NV_Ith_S (y_alpha, 2);
	const double H       = NV_Ith_S (y_alpha, 3);

	const double dV      = Bianchi_dVg_dphi (phi, c1->V0, c1->bV, c1->q);
	const double P       = Bianchi_P (alpha, a, phi, varphi, H, c1);
	const double D       = Bianchi_D (alpha, a, phi, varphi, H, c1);
	const double ddphi   = - (D * varphi + dV) / P;
	const double dH      = Bianchi_dH (alpha, a, phi, varphi, H, c1);

	NV_Ith_S (ydot_alpha, 0)   =  - a;
	NV_Ith_S (ydot_alpha, 1)   = - varphi / H;
	NV_Ith_S (ydot_alpha, 2)   = - ddphi / H;
	NV_Ith_S (ydot_alpha, 3)   = - dH / H;

	return 0;
}

/*******************************************************
*************** Integral: alpha ************************
******************* t(alpha) ***************************/

static int
Bianchi_t(realtype alpha, N_Vector y_alpha, N_Vector yQdot_alpha, void *f_data)
{
	const double H				= NV_Ith_S (y_alpha, 3);

	if(H<0)
	{
	NV_Ith_S (yQdot_alpha, 0) = 1.0 / fabs(H);
	}
	else
	{
		NV_Ith_S (yQdot_alpha, 0) =   - 1.0 / fabs(H);
	}

	return 0;	
}

/*******************************************************
****************** Root find: t ************************
********************************************************/

static int 
Bianchi_roots_t (realtype t, N_Vector y_t, realtype *gout, void *user_data)
{
	BianchiConsts *c1	= (BianchiConsts *) user_data;
	const double phi	= NV_Ith_S (y_t, 1);
	const double g      = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);
	
	
	gout[0] = NV_Ith_S (y_t, 3);
	gout[1] = NV_Ith_S (y_t, 2);
	gout[2] = g -1.0;
	
	return 0;
}

/*******************************************************
***************** Root find: alpha *********************
********************************************************/

static int 
Bianchi_roots_alpha (realtype alpha, N_Vector y_alpha, realtype *gout, void *user_data)
{
	BianchiConsts *c1		= (BianchiConsts *) user_data;
	const double phi	= NV_Ith_S (y_alpha, 1);
	const double g      = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);
	
	
	gout[0] = NV_Ith_S (y_alpha, 3);
	gout[1] = NV_Ith_S (y_alpha, 2);
	gout[2] = g -1.0;
	
	return 0;
}

/*******************************************************
******************** Run_evol: t ***********************
********************************************************/

SyState
run_evol_t(void *cvode,
			SyState init_t,
			BianchiConsts *c1,
			FILE *out[5],
			N_Vector y_t,
			N_Vector y0_t,
			const double sign_beta_p,
			const double sign_beta_m,
			const double eps,
			const double a_init,
			unsigned int dim,
			double reltol,
			double abstol,
			double alpha_stop
			//~ double **alpha_a,
			//~ double **a_powm3_a,
			//~ double **a_a,
			)
{
	printf("Entering run_evol_t\n");
	
	SyState init_alpha;
	int flag;
	double alpha, a, phi, varphi, H, t_p;
	
	int cont 		= 0;
	int cont1		= 0;

	double t_f		= 1.0e10;
	double t_i		= init_t.t_0;
	double a_i		= init_t.a;
	//~ double alpha_i	= init_t.alpha;
	//~ double chi_bV	= 1.0/(pow(c1->bV,(-c1->bV)/(c1->bV + 1)) + pow(c1->bV ,(1.0)/(c1->bV + 1)));
	double last_t	= t_i;
	
	NV_Ith_S (y0_t, 0) = init_t.a;
	NV_Ith_S (y0_t, 1) = init_t.phi;
	NV_Ith_S (y0_t, 2) = init_t.varphi;
	NV_Ith_S (y0_t, 3) = init_t.H;

/******************** Inicializing CVODE ***********************/

	flag = CVodeReInit (cvode, t_i, y0_t);
	CVODE_CHECK (&flag, "CVodeInit", 1, init_alpha);

	flag = CVodeSStolerances (cvode, reltol, abstol);
	CVODE_CHECK (&flag, "CVodeSStolerances", 1, init_alpha);

	flag = CVodeSetUserData (cvode, c1);
	CVODE_CHECK (&flag, "CVodeSetUserData", 1, init_alpha);

	flag = CVDense (cvode, dim);
	CVODE_CHECK (&flag, "CVDense", 1, init_alpha);

	flag = CVodeSetStopTime (cvode, t_f);
	CVODE_CHECK (&flag, "CVodeSetStopTime", 1, init_alpha);
	
	while (TRUE)
	{
		double t_step;

		flag = CVode (cvode, t_f, y_t, &t_step, CV_ONE_STEP /* CV_NORMAL */);
		CVODE_CHECK (&flag, "CVode", 1, init_alpha);
			
		double Hp, Hm;
						
		a      = NV_Ith_S (y_t, 0);
		phi    = NV_Ith_S (y_t, 1);
		varphi = NV_Ith_S (y_t, 2);
		H      = NV_Ith_S (y_t, 3);

		alpha = -log (a);

		Bianchi_H (t_step, a, phi, varphi, c1, &Hp, &Hm);
			
		//~ const double R			= Bianchi_R(a, H, c1);
		//~ const double B			= Bianchi_B(a, H, c1);
		const double g          = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);
		const double V          = Bianchi_Vg (phi, c1->V0, c1->bV, c1->q);
		const double rho        = Bianchi_rho (alpha, a, phi, varphi, H, c1);
		const double rho_sigma  = 0.5 * Bianchi_sigma2 (a, c1);
		const double p          = Bianchi_p (alpha, a, phi, varphi, H, c1);
		const double w          = p / rho ;
		const double H2         = H * H;
		const double Erro       = (H2 - (1.0/3.0) * rho - (1.0/6.0) * Bianchi_sigma2 (a, c1))/H2;
		const double rho_p		= rho + p;
		const double dH     	= Bianchi_dH (t_step, a, phi, varphi, H, c1);
		const double dda		= dH + H2;
		//~ const double dV      = Bianchi_dVg_dphi (phi, c1->V0, c1->bV, c1->q);
		//~ const double P       = Bianchi_P (alpha, a, phi, varphi, H, c1);
		//~ const double D       = Bianchi_D (alpha, a, phi, varphi, H, c1);
		//~ const double ddphi   = - (D * varphi + dV) / P;
			
		if (flag == CV_ROOT_RETURN)
		{
			int roots[2];
			
			flag = CVodeGetRootInfo (cvode, roots);
			CVODE_CHECK (&flag, "CVodeGetRootInfo", 1,  init_alpha);
		
			if ((roots[0])&&(dH > 0))
			{
				printf ("#Bounce found! % 20.15e\n", dH);
				
				if (out[8] != NULL)
				{
					fprintf (out[8], "tb = % 20.15e \n", 
					t_step
					);
				 }
				
			}
			
			if (roots[1])
			{
				printf ("#Turning point found % 20.15e\n", alpha);
				if (out[1] != NULL)
				{
					fprintf (out[8], "tp = % 20.15e\n", 
							t_step
							);
				 }
			}
			
			if (roots[2])
			{
				printf ("#Spurious roots found\n");
				
			}
					
		}
			
		if((fabs(g-1.0)<1.0e-7)&&(cont1==0))
		{
			/***activate this option to save de ghost condensate begining and ending times***/
			printf("# Ghost condesate begins % 20.15e", phi);
			//~ 
			if (out[0] != NULL)
			{
				fprintf (out[8], "tgh = % 22.15g", 
				t_step);
			 }
			 /********************/
			 
			 cont1 = cont1+ 1;
		}
			
		if (rho_p<0)
		{
			if (out[0] == NULL)
			{
				printf ("#B % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e\n",
				t_step, 
				a, 
				H,
				dda,
				Erro
				);
			}
		}
		else
		{
			if (out[0] == NULL)
			{
				printf ("#S % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e\n",
				t_step, 
				a, 
				H,
				dda,
				Erro
				);
			}
			
		}
		
		if (out[0] != NULL)
		{
			fprintf (out[0], "% 22.15e % 22.15e\n", 
			t_step, 
			H);
			
			fprintf (out[1], "% 22.15e % 22.15e\n", 
			t_step, 
			a);
			
			fprintf (out[2], "% 22.15e % 22.15e\n", 
			t_step, 
			rho_sigma);
			
			fprintf (out[3], "% 22.15e % 22.15e\n", 
			t_step, 
			rho);
			
			fprintf (out[4], "% 22.15e % 22.15e\n", 
			t_step, 
			V);
			
			fprintf (out[5], "% 22.15e % 22.15e\n", 
			t_step, 
			g);
			
			fprintf (out[6], "% 22.15e % 22.15e\n", 
			t_step, 
			w);
			
			fprintf (out[7], "% 22.15e % 22.15e\n", 
			t_step, 
			Erro);
	  }
		 
		/********** Activate this option to re-enter run_evol_alpha after ghost condensate *******/
				
			//~ if((rho_p > 0)&&(init_t.i == 0))
			//~ {
				//~ printf("Entrei na condição de parada rho+p>0\n");	
				//~ t_p = t_step;		
				//~ break;
			//~ }
		/****************************************************************************************/
		
		if (t_step == t_f)
		{
			printf("Insuficient t_f  in run_evol_t");
			t_p = t_step;
			break;
		}
			
		if (((fabs(t_step - last_t)/fabs(t_step))< 1.0e-15)&&(init_t.i == 1)&&(H>0.0))
		{
			printf("Saturated t.\n The end.\n");
			cont = cont+1;
			break;
		}
		else
		{
			last_t = t_step;
		}
			
		if((H>0)&&(a > 10* a_i))
		{
			printf("END: alpha = alpha_i, t = % 20.15e ", t_step);
			cont = cont + 1;
			break;
		}
			
		if((H<0)&&(a/a_i < 1.0e-3))
		{
			printf("Singularidade alpha = alpha_i, t = % 20.15e ", t_step);
			cont = cont+1;
			break;
		}
	
	}
	
	init_alpha.a        	= a;
	init_alpha.phi      	= phi;
	init_alpha.varphi   	= varphi;
	init_alpha.H        	= H;
	init_alpha.alpha_i   	= alpha;
	init_alpha.t_0			= t_p;
	init_alpha.cont			= cont;	
		
	return init_alpha;
}

/*******************************************************
******************** Run_evol: alpha *******************
********************************************************/

SyState
run_evol_alpha(void *cvode,
			SyState init_alpha,
			double alpha_i,
			BianchiConsts *c1,
			FILE *out[5],
			N_Vector y_alpha,
			N_Vector y0_alpha,
			N_Vector yQdot0_alpha,
			N_Vector yQdot_alpha,
			const double sign_beta_p,
			const double sign_beta_m,
			const double eps,
			const double a_init,
			unsigned int dim,
			double reltol,
			double abstol,
			double MaxStep,
			double alpha_stop
			)
{
	printf("Entering run_evol_alpha\n");
	
	SyState init_t;	
	int flag;
		
	double alpha_f, a, phi, varphi, H, t;
	
	int cont	= 0;
	int cont1	= 0;
	
	double last_alpha	= alpha_i;
	double a_i			= init_alpha.a;
	//~ double chi_bV		= 1.0/(pow(c1->bV,(-c1->bV)/(c1->bV + 1)) + pow(c1->bV,(1.0)/(c1->bV + 1)));
	
	NV_Ith_S (y0_alpha, 0)		= init_alpha.a;
	NV_Ith_S (y0_alpha, 1)		= init_alpha.phi;
	NV_Ith_S (y0_alpha, 2)		= init_alpha.varphi;
	NV_Ith_S (y0_alpha, 3)		= init_alpha.H;
	NV_Ith_S (yQdot0_alpha, 0)	= init_alpha.t_0;

/******************** Run_evol: alpha *******************/

	flag = CVodeReInit (cvode, alpha_i, y0_alpha);
	CVODE_CHECK (&flag, "CVodeInit", 1, init_t);
	//~ 
	flag = CVodeQuadReInit (cvode, yQdot_alpha);
	CVODE_CHECK (&flag, "CVodeQuadInit", 1, init_t);
	
	flag = CVodeSStolerances (cvode, reltol, abstol);
	CVODE_CHECK (&flag, "CVodeSStolerances", 1, init_t);

	flag = CVodeSetUserData (cvode, c1);
	CVODE_CHECK (&flag, "CVodeSetUserData", 1, init_t);

	flag = CVDense (cvode, dim);
	CVODE_CHECK (&flag, "CVDense", 1, init_t);

	//flag = CVodeSetStopTime (cvode, alpha_f);
	alpha_f = alpha_i + ((init_alpha.H > 0) ? - 100.0 : 100.0);
	flag = CVodeSetStopTime (cvode, alpha_f);
	CVODE_CHECK (&flag, "CVodeSetStopTime", 1, init_t);
	
	flag = CVodeSetMaxStep (cvode, MaxStep);
	CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, init_t);
	
	//~ flag = CVodeRootInit(cvode, 3, &Bianchi_roots_alpha); 
	//~ CVODE_CHECK (&flag, "CVodeSetStopTime", 1, init_t);

	printf("Goin from % 20.15e TO % 20.15e (H = % 20.15e)\n", alpha_i, alpha_f, init_alpha.H);
		
	while (TRUE)
	{	
		double alpha_step, Hp, Hm;
		
		flag = CVode (cvode, alpha_f, y_alpha, &alpha_step, CV_ONE_STEP /* CV_NORMAL */);
		CVODE_CHECK (&flag, "CVode", 1, init_t);
		
		flag = CVodeGetQuad (cvode, &alpha_step, yQdot_alpha);
		CVODE_CHECK (&flag, "CVodeGetQuad", 1, init_t);
 
		t		= NV_Ith_S (yQdot_alpha, 0);
		a		= NV_Ith_S (y_alpha, 0);
		phi		= NV_Ith_S (y_alpha, 1);
		varphi	= NV_Ith_S (y_alpha, 2);
		H		= NV_Ith_S (y_alpha, 3);

		
		Bianchi_H (alpha_step, a, phi, varphi, c1, &Hp, &Hm);

		const double rho        = Bianchi_rho (alpha_step, a, phi, varphi, H, c1);
		const double p          = Bianchi_p (alpha_step, a, phi, varphi, H, c1);
		const double w          = p / rho ;
		const double rho_sigma  = 0.5 * Bianchi_sigma2 (a, c1);
		const double g          = Bianchi_Vg (phi, c1->g0, c1->bg, c1->p);
		const double V          = Bianchi_Vg (phi, c1->V0, c1->bV, c1->q);
		const double H2         = H * H;
		const double Erro       = (H2 - (1.0/3.0) * rho - (1.0/6.0) * Bianchi_sigma2 (a, c1))/H2;
		const double rho_p		= rho + p;
		//~ const double R			= Bianchi_R(a, H, c1);
		//~ const double B			= Bianchi_B(a, H, c1);
		const double dH     	= Bianchi_dH (t, a, phi, varphi, H, c1);
		const double dda		= dH + H2;

		alpha_f = alpha_step;
			
	if (flag == CV_ROOT_RETURN)
	{
		int roots[2];
		flag = CVodeGetRootInfo (cvode, roots);
		CVODE_CHECK (&flag, "CVodeGetRootInfo", 1,  init_t);
	//~ 
		if (roots[0])
		{
			printf ("#Bounce found!\n");					
		}
		
		if (roots[1])
		{
			printf ("#Turning point found % 20.15e\n", alpha_step);
			if (out[1] != NULL)
			{
				fprintf (out[8], "tp = % 20.15e\n",
				t
				);
			 }
		}
		
		if (roots[2])
		{
			printf ("#Spurious roots found\n");
		}
				
	}
		
		if (out[0] == NULL)
			{
				printf ("#E % 22.15e % 22.15e % 22.15e % 22.15e % 22.15e\n",
				t, 
				a, 
				H,
				dda,
				Erro
				);
			}
		
		
			if (out[0] != NULL)
			{
				fprintf (out[0], "% 22.15e % 22.15e\n", 
			t, 
			H);
			
			fprintf (out[1], "% 22.15e % 22.15e\n", 
			t, 
			a);
			
			fprintf (out[2], "% 22.15e % 22.15e\n", 
			t, 
			rho_sigma);
			
			fprintf (out[3], "% 22.15e % 22.15e\n", 
			t, 
			rho);
			
			fprintf (out[4], "% 22.15e % 22.15e\n", 
			t, 
			V);
			
			fprintf (out[5], "% 22.15e % 22.15e\n", 
			t, 
			g);
			
			fprintf (out[6], "% 22.15e % 22.15e\n", 
			t, 
			w);
			
			fprintf (out[7], "% 22.15e % 22.15e\n", 
			t, 
			Erro);
			 }
		 
		 
		if((fabs(g-1.0)<1.0e-7)&&(cont1==0))
			{
				printf("Ghost condensate begins % 20.15e", phi);
				if (out[0] != NULL)
				{
					fprintf (out[8], "tgh = % 22.15g \n", 
					t    
					);
				 }
				 cont1 = cont1+ 1;
			}

		if(rho_p < 0)
		{
			init_t.alpha = alpha_step;
			printf("Condition rho+p < 0. End 'run_alpha'.\n");
			break;
		}
		
		if((H > 0)&&( a > 10 * init_alpha.a))
		{
			printf("Final scale factor equals to initial scale factor during expansion. End.");
			break;
		}				
	
		if ((H>0)&&( fabs(Erro)>1.0e-9 /*(fabs(alpha_step - last_alpha)/fabs(alpha_step))< 1.0e-10)*/)&&(rho_p > 0))
		{
			init_t.i		= 1;
			printf("Saturated alpha, breaking to go to run_evol_t, init_t.i %20.15d\n alpha %20.15e \n", init_t.i, last_alpha);
			break;
		}
		else
		{
			last_alpha = alpha_step;
			init_t.i		= 0;
		}

		//~ if((H>0)&&(fabs(alpha_step)>fabs(alpha_stop)))
		//~ {
			//~ printf("alpha = alpha_i, alpha = % 20.15e", alpha_step);
			//~ cont = cont+1;
			//~ break;
		//~ }
		
		if((H < 0)&&(a / a_i < 1.0e-3))
		{
			printf("Singularidade alpha =1e-3 alpha_i, t = % 20.15e ", t);
			cont = cont+1;
			break;
		}
		
		if(fabs(varphi) < 1.0e-8)
		{
			printf ("Turning point found % 20.15e\n", varphi);
					if (out[1] != NULL)
					{
						fprintf (out[8], "tp=% 20.15e", 
						t
						);
					 }
			
			}
	}
	
	init_t.a	        = a;
	init_t.alpha        = alpha_f;
	init_t.phi      	= phi;
	init_t.varphi   	= varphi;
	init_t.H        	= H;
	init_t.t_0			= t;
	init_t.cont			= cont;

	return init_t;
}


int 
main (int argc, char *argv[])
{
	printf("Entering MAIN\n");
	
	void *cvode_alpha = CVodeCreate (CV_BDF, CV_NEWTON);
	void *cvode_t = CVodeCreate (CV_BDF, CV_NEWTON);
	
	unsigned int dim         = 4;
	
	double       sign_beta_p = 0;
	double       sign_beta_m = 0;
	double       eps         = 0.0;
	
	double       alpha_i , alpha_f;
	double       t_i         = 0.0;
	double       abstol      = 0.0;
	double       reltol      = 1.0e-15;
	
	SyState init_alpha;
	SyState init_t;
	
	int flag;
	
	char format1[100];
	char format2[100];
	char format3[100];
	char format4[100]; 
	char format5[100];
	char format6[100];
	char format7[100];
	char format8[100];
	char format9[100];
	char format10[100];
	FILE *out[10] = {NULL, };

	BianchiConsts c1;

/******** Model parameters **********/

	c1.gamma    = 1.0e-3;
	c1.beta     = 5.0;
	c1.g0       = 1.1;
	c1.V0       = -1.0e-7;
	c1.bV       = 5.0;
	c1.bg       = 0.5;
	c1.p        = 0.01;
	c1.q        = 0.1;
	c1.sigma2_0 = 5.0e48;

/***********************************/

	switch (argc)
	{
		case 4:
		break;
		
		case 5:
		sprintf (format1, "%s-H-t.dat", argv[4]);
		sprintf (format2, "%s-a-t.dat", argv[4]);
		sprintf (format3, "%s-rho_sig-t.dat", argv[4]);
		sprintf (format4, "%s-rho_phi-t.dat", argv[4]);
		sprintf (format5, "%s-V-t.dat", argv[4]);
		sprintf (format6, "%s-g-t.dat", argv[4]);
		sprintf (format7, "%s-w-t.dat", argv[4]);
		sprintf (format8, "%s-E-t.dat", argv[4]);
		sprintf (format9, "%s-TPS.dat", argv[4]);
		sprintf (format10, "%s-PAR.dat", argv[4]);
		out[0] = fopen (format1, "w");
		out[1] = fopen (format2, "w");
		out[2] = fopen (format3, "w");
		out[3] = fopen (format4, "w");
		out[4] = fopen (format5, "w");
		out[5] = fopen (format6, "w");
		out[6] = fopen (format7, "w");
		out[7] = fopen (format8, "w");
		out[8] = fopen (format9, "w");
		out[9] = fopen (format10, "w");
 		break;
	
		default:
		fprintf (stderr, "usage: bianchi_phi eps sign_beta_p sign_beta_m OR bianchi_phi eps sign_beta_p sign_beta_m file_name to save the data.\n");
		exit (-1);
		break;
	}
	
	N_Vector y0_alpha		= N_VNew_Serial (dim);
	N_Vector yQdot0_alpha	= N_VNew_Serial (dim);
	N_Vector y_alpha		= N_VNew_Serial (dim);
	N_Vector yQdot_alpha	= N_VNew_Serial (dim);

	N_Vector y0_t = N_VNew_Serial (dim);
	N_Vector y_t  = N_VNew_Serial (dim);
	


			
/******* Initial conditions ******/
		double a_i = 1.0e10;	
		double phi_i		= -2.5;
		double varphi_i		= 8.0e-6;
/*************************************/	

	alpha_i =  - log(a_i);
	alpha_f	=  1000.0;
		
	init_alpha.alpha_i      = alpha_i;
	init_alpha.alpha_f      = alpha_f;
	init_alpha.t_0			= 0.0;
	
	init_alpha.a			= a_i;
	init_alpha.phi			= phi_i;
	init_alpha.varphi		= varphi_i;
	init_alpha.alpha		= alpha_i;

	eps         = atof (argv[1]);
	sign_beta_p = atof (argv[2]);
	sign_beta_m = atof (argv[3]);
		
	double Hp_i = 0.0;
	double Hm_i = 0.0;

	Bianchi_H (alpha_i, init_alpha.a, init_alpha.phi, init_alpha.varphi, &c1, &Hp_i, &Hm_i);
	
	double H_i;

	printf ("# Two possible constraints % 20.15e % 20.15e\n", Hp_i, Hm_i);
		
	/****************** Calculating H_i ********************/
		{
			H_i = Hp_i > Hm_i ? Hm_i : Hp_i;
			
			if (H_i > 0.0)
			{
				fprintf (stderr, "Bad initial conditions, both roots are positive.\n");
				exit (-1);
			}
	  
			if ((eps< 0.0) | ( eps>1.0))
			{
				fprintf (stderr, "0 < eps <=1");
				exit (-1);
			}

			if ( (fabs(sign_beta_p)!=1 ) | (fabs(sign_beta_m)!=1))
			{
				fprintf (stderr, "Signal atribution is: 1 for positive and -1 for negative");
				exit (-1);
			}       

		}
		
		init_alpha.H 	= H_i;
			
			 
/************** Initializing alpha *******************/
	
	flag = CVodeInit (cvode_alpha, &Bianchi_f_alpha, alpha_i, y0_alpha);
	CVODE_CHECK (&flag, "CVodeInit", 1, -1);
	
	flag = CVodeQuadInit (cvode_alpha, &Bianchi_t, yQdot0_alpha);
	CVODE_CHECK (&flag, "CVodeQuadInit", 1, -1);

/************** Initializing t ***********************/	

	flag = CVodeInit (cvode_t, &Bianchi_f_t, t_i, y0_t);
	CVODE_CHECK (&flag, "CVodeInit", 1, -1);
	
	flag = CVodeRootInit(cvode_t, 3, &Bianchi_roots_t); 
	CVODE_CHECK (&flag, "CVodeSetStopTime", 1, -1);
	
	flag = CVodeRootInit(cvode_alpha, 3, &Bianchi_roots_alpha); 
	CVODE_CHECK (&flag, "CVodeSetStopTime", 1, -1);
		
/******************************************************/		

	//~ double 	phi_lim, chi_bV;
			
	while(TRUE)
	{
		printf("Initization of 'run_alpha' \n t_i = % 20.15e \na_i = % 20.15e \nalpha_i = % 20.15e \nphi_i = % 20.15e \nvarphi_i = % 20.15e \nH_i = % 20.15e\n",
		init_alpha.t_0, init_alpha.a, init_alpha.alpha, init_alpha.phi, init_alpha.varphi, init_alpha.H);
	
		init_t = run_evol_alpha (cvode_alpha,
								init_alpha,
								init_alpha.alpha_i,
								&c1,
								out,
								y_alpha,
								y0_alpha,
								yQdot0_alpha,
								yQdot_alpha,
								sign_beta_p,
								sign_beta_m,
								eps,
								a_i,
								dim,
								reltol,
								abstol,
								0.0,
								alpha_i
								);
		
		if(init_t.cont==1)
		{
			break;
		}
		
		printf("Initization of 'run_t' \nt_i = % 20.15e \na_i = % 20.15e \nalpha_i = % 20.15e \nphi_i = % 20.15e \nvarphi_i = % 20.15e \nH_i = % 20.15e\n",
		init_t.t_0, init_t.a, init_t.alpha, init_t.phi, init_t.varphi, init_t.H);
		
		init_alpha = run_evol_t(
								cvode_t,
								init_t,
								&c1,
								out,
								y_t,
								y0_t,
								sign_beta_p,
								sign_beta_m,
								eps,
								 a_i,
								 dim,
								 reltol,
								 abstol,
								 alpha_i
								 );
		if(init_alpha.cont==1)
		{
			break;
		}
		
	}

/********************************************
************ SAVING PARAMETERS***************
********************************************/
	{
		//~ phi_lim = (sqrt((2.0)/(c1.q)) * log(1.0/c1.bV))/(c1.bV + 1);
		//~ chi_bV	= 1.0/(pow(c1.bV,(-c1.bV)/(c1.bV + 1)) + pow(c1.bV,(1.0)/(c1.bV + 1))); 

		if (out[0] != NULL)
		{
			fprintf (out[9], "phi_i =% 20.15g\n varphi_i=% 20.15g\n sigma_i^2 = % 20.15g\n H_i= % 20.15g\n g0 = % 20.15g\n V0 = % 20.15g\n p = % 20.15g\n q =  % 20.15g\n beta = % 20.15g\n gamma =  % 20.15g\n bV =  % 20.15g\n bg =  % 20.15g\n", 
			phi_i,
			varphi_i,
			c1.sigma2_0,
			H_i,
			c1.g0,
			c1.V0,
			c1.p,
			c1.q,
			c1.beta,
			c1.gamma,
			c1.bV,
			c1.bg			
			);
		}
	}
  
	{
		int i;
		
		for (i=0;i<10;i++)
		{
			if (out[i] != NULL)
			fclose (out[i]);
		}
	}
  
  N_VDestroy (y0_t);
  N_VDestroy (y_t);
  N_VDestroy (y0_alpha);
  N_VDestroy (y_alpha);
  CVodeFree (&cvode_alpha);
  CVodeFree (&cvode_t);

  return (0);
}

int 
cvode_check_flag (void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  switch (opt)
  {
    case 0:
    {
      if (flagvalue == NULL)
      {
        fprintf (stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return FALSE;
      }
      break;
    }
    case 1:
    {
      errflag = (int *) flagvalue;
      if (*errflag < 0)
      {
        fprintf (stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
        return FALSE;
      }
      break;
    }
    case 2:
    {
      if (flagvalue == NULL)
      {
        fprintf (stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return FALSE;
      }
      break;
    }
    default:
    {
      fprintf (stderr, "Error: should not be reached!\n");
      exit (-1);
    }
  }
  return TRUE;
}

double 
sqrt1px_m1 (const double x)
{
  double binfact = 1.0;
  double res = 0;
  double xn = 1;
  int n = 0;

  if (fabs(x) > 1e-1)
    return sqrt (1.0 + x) - 1.0;

  while (TRUE)
  {
    binfact *= 3.0 / (2.0 * (1.0 + n++)) - 1.0;
    xn *= x;
    res += binfact * xn;

    if (fabs (binfact * xn / res) < DBL_EPSILON)
      break;
  }

  return res;
}
