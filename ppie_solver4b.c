#include "functions_method7.h"
#include <errno.h> 
#include <string.h>


struct moments
{
	double favg; // <f>
	double f2avg; // <f^2>
	double sfavg; // <sf>
	double fefiavg; // <fe fi>
	double chi2; // \chi^2
};

int integraleqn_iterator_parallel(int iteration, int seed, double *D, void *femoments, void *fimoments, double *offsets, void * params)
{

	struct pdf_params * modelparamstemp = (struct pdf_params *) params;
	struct moments * fe_moments = (struct moments *) femoments;
	struct moments * fi_moments = (struct moments *) fimoments;

	struct pdf_params modelparams;

	modelparams.stimavg=0, modelparams.sigstim = modelparamstemp->sigstim; 
	modelparams.avgin=0, modelparams.sigin = modelparamstemp->sigin, modelparams.vin = modelparamstemp->vin; 
	modelparams.avgout=0, modelparams.sigout = modelparamstemp->sigout, modelparams.vout = modelparamstemp->vout;
	modelparams.Poissk=modelparamstemp->Poissk;
	
	// Make more convenient names.
	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	double Poissk = modelparams.Poissk;
	//double z = modelparams->z; not needed right now
	


	// Define flag to indicate if iteration has converged.
	int doneflag = 1; // Flag to indicate if iteration is complete. If any of the bounds are not satisfied, the flag is set to 0.

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Define model parameter groups, vectors, etc.
	//
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double A = sigstim*sigstim/(sigstim*sigstim+sigin*sigin); // Combination of parameters that appears in int. eqn.

	double scale = sqrt(2*(sigstim*sigstim + sigin*sigin));
	double veff = (sigstim*sigstim+sigin*sigin*vin)/(sigstim*sigstim+sigin*sigin);
	double Lambda2_diff = sigout*sigout*(1-vout);
	double Lambda2_sum = sigout*sigout*(1+vout);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Import data defining functions fe and fi at several points and construct interpolating function.     
	//
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// The functions to be interpolated are not necessarily within the bounds we wish to specify (0 < f < 1). We
	// will interpolate the entire function we are given and use our interpolant to find where it crosses our 
	// desired threholds.

	// Data filename set-up

	FILE *feinput;
	FILE *fiinput;
	FILE *feoutput;
	FILE *fioutput;
	FILE *Doutput;

	char feinname[200], fiinname[200]; // input file name arrays
	char feoutname[200], fioutname[200]; // output file name arrays
	char Doutname[200];	

	sprintf(feinname,"%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%i%s%i%s","../data/fe_function_sigstim",sigstim,"_sigin",sigin,"_sigout",sigout,"_vin",vin,"_vout",vout,"_kappa",Poissk,"_veff",veff,"_seed",seed,"_iter",iteration,".txt");
	sprintf(fiinname,"%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%i%s%i%s","../data/fi_function_sigstim",sigstim,"_sigin",sigin,"_sigout",sigout,"_vin",vin,"_vout",vout,"_kappa",Poissk,"_veff",veff,"_seed",seed,"_iter",iteration,".txt");

	sprintf(feoutname,"%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%i%s%i%s","../data/fe_function_sigstim",sigstim,"_sigin",sigin,"_sigout",sigout,"_vin",vin,"_vout",vout,"_kappa",Poissk,"_veff",veff,"_seed",seed,"_iter",iteration+1,".txt");
	sprintf(fioutname,"%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%i%s%i%s","../data/fi_function_sigstim",sigstim,"_sigin",sigin,"_sigout",sigout,"_vin",vin,"_vout",vout,"_kappa",Poissk,"_veff",veff,"_seed",seed,"_iter",iteration+1,".txt");
	sprintf(Doutname,"%s","../data/Dcoeffs.txt");


	// Set variables to record data
	int fecount=0, ficount=0; // Counts how many rows are in each data file.
	float fedummy1=0, fedummy2=0; // Dummy variables.
	float fidummy1=0, fidummy2=0;

	feinput = fopen(feinname,"r");
	fiinput = fopen(fiinname,"r");
	
	// Read in from files (first run) to count number of points
	while(fscanf(feinput,"%f\t %f\n",&fedummy1,&fedummy2) == 2)
		fecount++;

	fclose(feinput);

	while(fscanf(fiinput,"%f\t %f\n",&fidummy1,&fidummy2) == 2)
		ficount++;

	fclose(fiinput);

	// Dynamically allocate memory to store data

	double *ze_data = (double *) malloc(fecount*sizeof(double)); // z-coordinate data for fe(z)
	double *zi_data = (double *) malloc(ficount*sizeof(double)); // z-coordinate data for fi(z)
	double *fe_data = (double *) malloc(fecount*sizeof(double)); // values of fe(z)
	double *fi_data = (double *) malloc(ficount*sizeof(double)); // values of fi(z)

	// Read in from files again and store in data arrays.

	feinput = fopen(feinname,"r");
	fiinput = fopen(fiinname,"r");

	int j = 0; // temp counter j.
	while(fscanf(feinput,"%f\t %f\n",&fedummy1,&fedummy2) == 2)
	{
		ze_data[j] = fedummy1;
		fe_data[j] = fedummy2;
		j++;
	}

	fclose(feinput);

	j = 0; // reset counter j.
	while(fscanf(fiinput,"%f\t %f\n",&fidummy1,&fidummy2) == 2)
	{
		zi_data[j] = fidummy1;
		fi_data[j] = fidummy2;
		j++;
	}

	fclose(fiinput);

	// Begin interpolation. Interpolate with cubic splines (in case there's any curvature).

	gsl_interp_accel *fe_acc = gsl_interp_accel_alloc ();
	gsl_spline *fe_spline = gsl_spline_alloc (gsl_interp_cspline, fecount);

	gsl_interp_accel *fi_acc = gsl_interp_accel_alloc ();
	gsl_spline *fi_spline = gsl_spline_alloc (gsl_interp_cspline, ficount);

	gsl_spline_init (fe_spline, ze_data, fe_data, fecount);

	gsl_spline_init (fi_spline, zi_data, fi_data, ficount);

	// To evaluate spline at a given z, use fe_eval = gsl_spline_eval(fe_spline, z, fe_acc). Similarly,
	// fi_eval = gsl_spline_eval(fi_spline, z, fi_acc).


	// Needs 'roots' for which the nonlinearities cross the constraint boundaries at 0 and 1.
	
	// The 'roots' are actually bounds on the roots for given decoder coefficients. 
	// The form depends on the sign of the decoder coefficients. We select the appropriate 
	// root bounds below.
	
	double ze_0 = (D[1]*Poissk/2.0 + D[0] + D[2]*heaviside(-D[1]*D[2]))/A;
	double ze_1 = (D[1]*(Poissk/2.0+1.0) + D[0] + D[2]*heaviside(D[1]*D[2]))/A;
	
	double zi_0 = -(D[2]*Poissk/2.0 + D[0] + D[1]*heaviside(-D[1]*D[2]))/A;
	double zi_1 = -(D[2]*(Poissk/2.0+1.0) + D[0] + D[1]*heaviside(D[1]*D[2]))/A;
	
//	printf("ze_0 = %g, ze_1 = %g, zi_0 = %g, zi_1 = %g\n",ze_0,ze_1,zi_0,zi_1);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Find offsets: values of z for which f(z) = 0.5.
	//
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////


	struct rootfindingparams fe_splineparams; // Define struct to hold spline data and desired root

	fe_splineparams.splineptr = fe_spline;
	fe_splineparams.accptr = fe_acc;
	fe_splineparams.root = 0.5;

	struct rootfindingparams fi_splineparams; // Define struct to hold spline data and desired root

	fi_splineparams.splineptr = fi_spline;
	fi_splineparams.accptr = fi_acc;
	fi_splineparams.root = 0.5;

	int introot=0; // for loop index
	double offset_e=0; // Define variables to hold results of root locations.
		
	int status_fe;
	int iter_fe = 0, max_iter_fe = 100;
	const gsl_root_fsolver_type *T_fe;
	gsl_root_fsolver *s_fe;
	double r_fe = 0;
	double x_lo_fe = ze_data[0], x_hi_fe = ze_data[fecount-1];
	gsl_function F_fe;

	F_fe.function = &interpolantminusxval; 
	F_fe.params = &fe_splineparams;

	T_fe = gsl_root_fsolver_brent;
	s_fe = gsl_root_fsolver_alloc (T_fe);
	gsl_root_fsolver_set (s_fe, &F_fe, x_lo_fe, x_hi_fe);

	do
	{
		iter_fe++;
		status_fe = gsl_root_fsolver_iterate (s_fe);
		r_fe = gsl_root_fsolver_root (s_fe);
		x_lo_fe = gsl_root_fsolver_x_lower (s_fe);
		x_hi_fe = gsl_root_fsolver_x_upper (s_fe);
		status_fe = gsl_root_test_interval (x_lo_fe, x_hi_fe,0.0001, 0.0001);
	}
	while (status_fe == GSL_CONTINUE && iter_fe < max_iter_fe);

	if(iter_fe == max_iter_fe)
		printf("WARNING! Root finder may not have converged!\n");
	else
	{
		//printf("fe root = %g\n",r_fe);
	}

	offset_e = r_fe; // Record fe(z) = 1 root.

	gsl_root_fsolver_free (s_fe);


	// Do for fi now.

	double offset_i=0; // Define variables to hold results of root locations.
	
	fi_splineparams.root = (double) introot;
	int status_fi;
	int iter_fi = 0, max_iter_fi = 100;
	const gsl_root_fsolver_type *T_fi;
	gsl_root_fsolver *s_fi;
	double r_fi = 0;
	double x_lo_fi = zi_data[0], x_hi_fi = zi_data[ficount-1];
	gsl_function F_fi;

	F_fi.function = &interpolantminusxval; 
	F_fi.params = &fi_splineparams;

	T_fi = gsl_root_fsolver_brent;
	s_fi = gsl_root_fsolver_alloc (T_fi);
	gsl_root_fsolver_set (s_fi, &F_fi, x_lo_fi, x_hi_fi);

	do
	{
		iter_fi++;
		status_fi = gsl_root_fsolver_iterate (s_fi);
		r_fi = gsl_root_fsolver_root (s_fi);
		x_lo_fi = gsl_root_fsolver_x_lower (s_fi);
		x_hi_fi = gsl_root_fsolver_x_upper (s_fi);
		status_fi = gsl_root_test_interval (x_lo_fi, x_hi_fi,0.0001, 0.0001);
			//printf("iter_fi = %i, x_lo_fi=%g, x_hi_fi=%g, r_fi = %g\n",iter_fi,x_lo_fi,x_hi_fi,r_fi);
	}
	while (status_fi == GSL_CONTINUE && iter_fi < max_iter_fi);

	if(iter_fi == max_iter_fi)
		printf("WARNING! Root finder may not have converged! (root = %g)\n",r_fi);
	else
	{
		//printf("fi root = %g\n",r_fi);
	}

	offset_i = r_fi; // Record fi(z) = 1 root.

	gsl_root_fsolver_free (s_fi);  

	*(offsets+0) = offset_e;
	*(offsets+1) = offset_i;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// We need to compute several integrals over these interpolated functions. We do that below.
	//
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Set up Gauss-Legendre integration table, workspaces, and pointers to be passed to integration functions.

	// We'll use separate tables for the E and I functions, though for our cases of interest a single table is
	// probably sufficient.
	
	int numpts = 3; // Want npts*fecount (npts*ficount) points in our integration tables to get sufficient resolution

	gsl_integration_glfixed_table * t_fe = gsl_integration_glfixed_table_alloc (numpts*fecount);
	gsl_integration_glfixed_table * t_fi = gsl_integration_glfixed_table_alloc (numpts*ficount);

	// We will pass these tables, along with the spline pointers and arrays containing the x-data for the 
	// functions, in a struct.

	struct interpolantintegrationparams intparams_fe;
	intparams_fe.GLtable = t_fe;
	intparams_fe.splineptr = fe_spline;
	intparams_fe.accptr = fe_acc;
	intparams_fe.numintpts = numpts;
	intparams_fe.zeroroot = ze_0;
	intparams_fe.unityroot = ze_1;
	intparams_fe.xcount = fecount;
	intparams_fe.modelparams = modelparams;

	struct interpolantintegrationparams intparams_fi;
	intparams_fi.GLtable = t_fi;
	intparams_fi.splineptr = fi_spline;
	intparams_fi.accptr = fi_acc;
	intparams_fi.numintpts = numpts;
	intparams_fi.zeroroot = zi_0;
	intparams_fi.unityroot = zi_1;
	intparams_fi.xcount = ficount;
	intparams_fi.modelparams = modelparams;

	struct fefi_doubleintegralparams doubleparams;
	doubleparams.GLtable_fe = t_fe;
	doubleparams.splineptr_fe = fe_spline;
	doubleparams.accptr_fe = fe_acc;
	doubleparams.zeroroot_fe = ze_0;
	doubleparams.unityroot_fe = ze_1;
	doubleparams.fe_count = fecount;
	doubleparams.GLtable_fi = t_fi;
	doubleparams.splineptr_fi = fi_spline;
	doubleparams.accptr_fi = fi_acc;
	doubleparams.zeroroot_fi = zi_0;
	doubleparams.unityroot_fi = zi_1;
	doubleparams.fi_count = ficount;
	doubleparams.numintpts = numpts;
	doubleparams.modelparams = modelparams;

	//printf("&intparams_fe = %p, &intparams_fi = %p\n",&intparams_fe,&intparams_fi);
	
	// Integral < fe(s + \eta_e) >
	
	double fe_sum = fe_avg(&intparams_fe);
	
	// Integral < fi(-s - \eta_e) >
	
	double fi_sum = fi_avg(&intparams_fi);


	// Integral: < s fe(s+\eta_e) > = \sigma_s^2/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{3/2} \times
	// \int_{-\infty}^\infty dx x \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fe(x)

	double sfe_sum = sfe_avg(&intparams_fe);

	// Integral: < s fi(-s-\eta_i) > = -\sigma_s^2/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{3/2} \times
	// \int_{-\infty}^\infty dx x \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fi(x)

	double sfi_sum = sfi_avg(&intparams_fi);

	// Integral: < f^2e(s+\eta_e) > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fe(x)^2

	double fe2_sum = fe2_avg(&intparams_fe);

	// Integral: < f^2i(-s+\eta_i) > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fi(x)^2

	double fi2_sum = fi2_avg(&intparams_fi);

	// Double integral <f_E(s+\eta_e) f_I(-s-\eta_I) >

	double fefi_sum = fegfi_montecarlo_integral(&doubleparams);

	// Set values in struct to pass back to min loop (chi^2 recorded next section)
	
	fe_moments->favg = fe_sum, fe_moments->f2avg = fe2_sum, fe_moments->sfavg = sfe_sum;
	fi_moments->favg = fi_sum, fi_moments->f2avg = fi2_sum, fi_moments->sfavg = sfi_sum;
	fe_moments->fefiavg = fefi_sum, fi_moments->fefiavg = fefi_sum; 

	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////
	// 
	// Compute the linear estimator coefficients De and Di.
	//
	///////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////

	// De = [(<zeta_i^2>+<fi^2>)*<s*fe> - (<zeta_e*zeta_i>+<fi*fe>)*<s*fi>]/[(<zeta_e^2>+<f_e^2>)(<zeta_i^2>+<f_i^2>)-(<zeta_e*zeta_i>+<fi*fe>)^2]

	// Di = [(<zeta_e^2>+<fe^2>)*<s*fi> - (<zeta_e*zeta_i>+<fi*fe>)*<s*fe>]/[(<zeta_e^2>+<f_e^2>)(<zeta_i^2>+<f_i^2>)-(<zeta_e*zeta_i>+<fi*fe>)^2]

	double covarout = sigout*sigout*vout;
	
	double var_fe = fe2_sum - fe_sum*fe_sum;
	double var_fi = fi2_sum - fi_sum*fi_sum;
	double covar_fefi = fefi_sum - fe_sum*fi_sum;

	double D_den = (sigout*sigout + var_fe + Poissk*fe_sum)*(sigout*sigout + var_fi + Poissk*fi_sum)-(covarout+covar_fefi)*(covarout+covar_fefi);

	double De_num = (sigout*sigout + var_fi + Poissk*fi_sum)*sfe_sum - (covarout + covar_fefi)*sfi_sum;

	double Di_num = (sigout*sigout + var_fe + Poissk*fe_sum )*sfi_sum - (covarout + covar_fefi)*sfe_sum;

	double Deold = D[1]; // Store previous values of De and Di so that we can check convergence.
	double Diold = D[2];

	double D0 = -(D[1]*fe_sum + D[2]*fi_sum); // Compute D0 (using values from step n, for consistency).
	double De = De_num/D_den;
	double Di = Di_num/D_den;

	// Updated root bounds are need FOR THE RANGE OF Z TO PRINT ONLY
	
	double ze_0_new = (De*Poissk/2.0 + D0 + Di*heaviside(-De*Di))/A;
	double ze_1_new = (De*(Poissk/2.0+1.0) + D0 + Di*heaviside(De*Di))/A;
	
	double zi_0_new = -(Di*Poissk/2.0 + D0 + De*heaviside(-De*Di))/A;
	double zi_1_new = -(Di*(Poissk/2.0+1.0) + D0 + De*heaviside(De*Di))/A;
	

	// Write new De, Di to D pointer.

	*(D+0) = D0;
	*(D+1) = De;
	*(D+2) = Di; 
 
	// Print D coeffs to file
	Doutput = fopen(Doutname,"a");
	fprintf(Doutput,"%i\t %g\t %g\t %g\n",iteration,D0,De,Di);
	fclose(Doutput);

	printf("D0 = %g, De = %g, Di = %g\n",D0,De,Di);

	// First convergence check

	double tol = 1e-4; // Set the tolerance for our solution to be considered converged. 
			    //i.e., want |f^(n+1)(z) - f^(n)(z)|/sigout < tol, 
		            //and |D^(n)-D^(n-1)| < tol.

	printf("|De - De| = %g, |Di - Di| = %g\n",fabs(De-Deold),fabs(Di-Diold));
	if(fabs(De-Deold) > tol || fabs(Di-Diold) > tol)
		doneflag = 0; 
		
	// Print value of mean square error \chi^2 = \sigma_s^2 + 2*D*b-4*D*<sf> - 2*D^2*(\kappa*<f> + <f^2>) - 2*D*<f_+ f_->
	//double chi2 = sigstim*sigstim - 2.0*(De*sfe_sum + Di*sfi_sum) + De*De*(sigout*sigout+ Poissk*fe_sum + var_fe) + Di*Di*(sigout*sigout+ Poissk*fi_sum + var_fi) +2.0*De*Di*(covarout + covar_fefi);

	double chi2 = sigstim*sigstim - De*sfe_sum - Di*sfi_sum; //chi^2 for optimal soln. Maybe not be accurate for intermediate iterations.

   // Print to output structs
   fe_moments->chi2 = chi2;
   fi_moments->chi2 = chi2; //redundant, but might as well...
   
//   printf("%g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n",Deold,Diold,fe_moments->favg,fi_moments->favg,fe_moments->f2avg,fi_moments->f2avg,fe_moments->sfavg, fi_moments->sfavg, fe_moments->fefiavg,fe_moments->chi2); //print moments & D's.
//   fprintf(stdout,"%g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n",Deold,Diold,fe_sum,fi_sum,fe2_sum,fi2_sum,sfe_sum, sfi_sum, fefi_sum, chi2); //print moments & D's.	
//   fprintf(stdout,"De=%g, Di=%g, <fe> = %g, <fi> = %g, <fe2> = %g, <fi2> = %g, <sfe> = %g, <sfi> = %g, <fefi> = %g, chi2 = %g\n",Deold,Diold,fe_sum,fi_sum,fe2_sum,fi2_sum,sfe_sum, sfi_sum, fefi_sum, chi2); //print moments & D's.	
	
//	fprintf(stdout,"test = %g\n",fe_constrained(-0.02,&intparams_fe));
	
//	fprintf(stdout,"favg test = %g\n",fe_avg(&intparams_fe));

	// Open output file names

	feoutput = fopen(feoutname,"w");
	fioutput = fopen(fioutname,"w");

	double feK1_sum=0, fiK1_sum=0;
	double fe = 0, fi = 0;
	double z=0;
	//double z0 = -3.0*maxval(fabs(De),fabs(Di))*(1+sigin*sigin/sigstim/sigstim)*(1+Poissk/2+ fabs(D0/minval(fabs(De),fabs(Di))));
	double zstep = 0.01; 
	double z0 = 3.0*maxval(maxval(fabs(ze_0_new),fabs(zi_0_new)),maxval(fabs(ze_1_new),fabs(zi_1_new)));
	double zmax = ceil(z0/zstep)*zstep; 
	double feold = 0, fiold = 0; // Variables to store old values of fe and fi


	for(z=-zmax;z<= zmax; z += zstep)
	{

		modelparams.z = z; // need to update value of z passed to upperintegratedK1 
		intparams_fe.modelparams = modelparams; // and and feK1_sum.
		intparams_fi.modelparams = modelparams; // and and feK1_sum.
		
		// Use old z bounds here since it refers to old guesses for f and D.
		
//		printf("checkpoint 1...\n");
		
		// compute \int_{ze_0}^{ze_1} dx~f_E(x)K1(-x,z)
		    feold = fe_constrained(z,&intparams_fe);
		    
//		printf("checkpoint 2...\n");
		    
			feK1_sum = feK1integral(&intparams_fe); 
			
//		printf("checkpoint 3...\n");

		
		// compute \int_{ze_0}^{ze_1} dx~f_E(x)K1(-x,z)
		    fiold = fi_constrained(z,&intparams_fi);
			fiK1_sum = fiK1integral(&intparams_fi); 

		//printf("fiold = %g\n",fiold);
		
		fe = (A*z - Di*fiK1_sum - D0)/De - Poissk/2;
		fi = -(A*z + De*feK1_sum + D0)/Di - Poissk/2;

		fprintf(feoutput,"%.20g\t %.20g\n",z,fe);
		fprintf(fioutput,"%.20g\t %.20g\n",z,fi);

		// Check for convergence
		if(doneflag == 1)
		{
			if(fabs(rangethreshold(fe,1.0) - feold) > tol || fabs(rangethreshold(fi,1.0)-fiold) > tol)
			{
				doneflag = 0;
				printf("%g\t, |fe-fe| = %g\t, |fi -fi| = %g\n",z,fabs(rangethreshold(fe,1.0) - feold),fabs(rangethreshold(fi,1.0)-fiold));
			}
		}

	}

	// Free interpolant constructs

	gsl_spline_free (fe_spline);
	gsl_interp_accel_free (fe_acc);
	gsl_spline_free (fi_spline);
	gsl_interp_accel_free (fi_acc);

	// Free integration tables

	gsl_integration_glfixed_table_free (t_fe);
	gsl_integration_glfixed_table_free (t_fi);

	// Close output files

	fclose(feoutput);
	fclose(fioutput);

	return doneflag;
}

int main(int argc, char *argv[])
{

	// Convert input values to convenient names.
	
	// Write input parameters to user-friendly names
	double sigstim = atof(argv[1]);
	double sigin = atof(argv[2]);
	double sigout = atof(argv[3]);
	double vin = atof(argv[4]);
	double vout = atof(argv[5]);
	double Poissk = atof(argv[6]); 
	int seed = atoi(argv[7]);
	int OOvOF = atoi(argv[8]); //ON-ON = 0, ON-OFF = 1
	
	if(OOvOF != 0 && OOvOF != 1)
		printf("Warning! Last argument must be 0 (ON-ON) or 1 (ON-OFF)\n");

	int iteration = 0; //which iteration of the function solution we are on. Will read current iteration in from file
	int max_iterations = 300; // Maximum number of iterations to do.

	struct pdf_params modelparams; // Define struct to contain parameters.

	// Set values
	modelparams.sigstim = sigstim;
	modelparams.sigin = sigin, modelparams.vin = vin; 
	modelparams.sigout = sigout, modelparams.vout = vout;
	modelparams.Poissk = Poissk;
	modelparams.z = 0; // This parameter is a variable to be used later; set to zero for now.
	
	double veff = (sigstim*sigstim + sigin*sigin*vin)/(sigstim*sigstim+sigin*sigin);

	// Generate and print first guesses, f_E(z) = A*z, f_I(z) = -A*z.

	FILE *feoutput;
	FILE *fioutput;
	FILE *iterationcounter;
	FILE *momentsoutput;

	char feoutname[200], fioutname[200], momentsoutname[200], itcounter[300]; // output file name arrays	

	sprintf(feoutname,"%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%i%s%i%s","../data/fe_function_sigstim",sigstim,"_sigin",sigin,"_sigout",sigout,"_vin",vin,"_vout",vout,"_kappa",Poissk,"_veff",veff,"_seed",seed,"_iter",0,".txt");
	sprintf(fioutname,"%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%i%s%i%s","../data/fi_function_sigstim",sigstim,"_sigin",sigin,"_sigout",sigout,"_vin",vin,"_vout",vout,"_kappa",Poissk,"_veff",veff,"_seed",seed,"_iter",0,".txt");
	sprintf(itcounter,"%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%i%s","../data/iterationcounter_sigstim",sigstim,"_sigin",sigin,"_sigout",sigout,"_vin",vin,"_vout",vout,"_kappa",Poissk,"_veff",veff,"_seed",seed,".txt");


	feoutput = fopen(feoutname,"w+");
	fioutput = fopen(fioutname,"w+");
	
	if(fioutput == NULL)
	{
		printf("fioutput = NULL\n");
	}
		
		

	double fe = 0, fi = 0;
	double z=0; 
	int doneflag = 0; //Define variable for flag which determines when iteration is finished.
	double * D;
	double * offsets; 
	D = malloc(3*sizeof(double));
	offsets = malloc(2*sizeof(double));
	D[0] = 0.0, D[1] = 0.95, D[2] = 1.0; // Define initial vector defining the linear coefficients D0 = D[0], De = D[1] and Di = D[2].
 	D[0] = -(D[1]+D[2])/2; // Initial guess.

	int itfileexists = 1;
	iterationcounter = fopen(itcounter,"r");		
	if(iterationcounter == NULL)
	{
		//perror(itcounter);
		itfileexists = 0;
	}
	
//	printf("itfileexists = %i\n",itfileexists);
	if(itfileexists == 1)
	{
		int itdummy=0;
		float D1dummy=0, D2dummy=0;
		//while(fscanf(iterationcounter,"%i\n",&itdummy) != EOF)
		//{
			fscanf(iterationcounter,"%i\t %f\t %f\n",&itdummy,&D1dummy,&D2dummy);
			iteration = itdummy;
			D[1] = D1dummy;
			D[2] = D2dummy;
			D[0] = -(D[1]+D[2])/2; // Initial guess; really -(D1*<f1> + D2*<f2>)
		//	printf("itdummy = %i\n",itdummy);
		//}
//		if(iteration > 0)
//			iteration = iteration - 1; // Set iteration back 1 in case most recent iteration did not finish
//		else
//			iteration = 0; // Make sure iteration is zero otherwise.
	}
	
	fclose(iterationcounter);
 	
 	 	// Set up structs to hold moments
 	
 	struct moments *fe_moments = malloc(sizeof(struct moments));
 	fe_moments->favg = 0, fe_moments->f2avg = 0, fe_moments->sfavg = 0, fe_moments->fefiavg = 0, fe_moments->chi2=0;
	
	struct moments *fi_moments = malloc(sizeof(struct moments));
	fi_moments->favg = 0, fi_moments->f2avg = 0, fi_moments->sfavg = 0, fi_moments->fefiavg = 0, fi_moments->chi2=0;
 	
	double svar = sigstim*sigstim, invar = sigin*sigin;
	double A = sigstim*sigstim/(sigstim*sigstim+sigin*sigin); // Combination of parameters that appears in int. eqn.
	double scale = sqrt(2*invar*(1-vin)*(2*svar+invar*(1+vin))/(svar+invar));
 		
	// Set up RNG to generate initial function guess.
	const gsl_rng_type *T;
  	gsl_rng *r;
    gsl_rng_env_setup();

  	T = gsl_rng_default;
  	r = gsl_rng_alloc (T);

    gsl_rng_set(r,seed);
  	
  	double ue=0, ui=0; // random vars.
  	double De_r=1, Di_r=1;
 		
 		
 	if(iteration == 0)	
 	{
 	
 	  	//u = gsl_rng_uniform (r);
  	
		De_r = 1 + gsl_ran_gaussian(r,0.5); //Gaussian rv randomly centered at +1.
		Di_r = 1 + gsl_ran_gaussian(r,0.5); 
		if(OOvOF == 1)
			Di_r = -Di_r; // OFF cell guess for neuron 2.

		D[1] = De_r; 
		D[2] = Di_r;

		D[0] = -(D[1]+D[2]);

 		double ze_0 = (D[1]*Poissk/2.0 + (fabs(D[1]) + 2.0*fabs(D[2]))*(2.0*heaviside(-D[1])-1.0))/A;
		double ze_1 = (D[1]*(Poissk/2.0+1.0) + (fabs(D[1]) + 2.0*fabs(D[2]))*(2.0*heaviside(D[1])-1.0))/A;
	
		double zi_0 = (-D[2]*Poissk/2.0 + (fabs(D[2]) + 2.0*fabs(D[1]))*(2.0*heaviside(D[2])-1.0))/A;
		double zi_1 = (-D[2]*(Poissk/2.0+1.0) + (fabs(D[2]) + 2.0*fabs(D[1]))*(2.0*heaviside(-D[2])-1.0))/A;
	
		printf("ze_0 = %g, ze_1 = %g, zi_0 = %g, zi_1 = %g\n",ze_0,ze_1,zi_0,zi_1);
 	
 	
 		//double z0 = -3.0*maxval(fabs(D[2]),fabs(D[1]))*(1+sigin*sigin/sigstim/sigstim)*(1+Poissk/2 + fabs(D[0]/minval(fabs(D[2]),fabs(D[1]))));
		double zstep = 0.01; 
		double z0 = 3.0*maxval(maxval(fabs(ze_0),fabs(zi_0)),maxval(fabs(ze_1),fabs(zi_1)));
		double zmax = ceil(z0/zstep)*zstep; 
	
 	

		for(z=-zmax;z<= zmax; z += zstep)
		{

			ue = (2.0*gsl_rng_uniform (r)-1.0)*(fabs(De_r) + 2.0*fabs(Di_r));
			ui = (2.0*gsl_rng_uniform (r)-1.0)*(fabs(Di_r) + 2.0*fabs(De_r));
		
			fe = (z + ue)/De_r - Poissk/2;
			fi = (-z + ui)/Di_r - Poissk/2;


			fprintf(feoutput,"%.20g\t %.20g\n",z,fe);
			fprintf(fioutput,"%.20g\t %.20g\n",z,fi);
		}

		fclose(feoutput);
		fclose(fioutput);
		
	}

	// Now, generate updated guess for De, Di, f_E(z) and f_I(z) using a Picard Iteration-like method.

	do
	{
		doneflag = integraleqn_iterator_parallel(iteration,seed,D,fe_moments,fi_moments,offsets,&modelparams);;
		iteration++;
		
		// Print current iteration state to file along with D1 and D2
		iterationcounter = fopen(itcounter,"w");
		fprintf(iterationcounter,"%i\t %g\t %g\n",iteration,D[1],D[2]);
		fclose(iterationcounter);

		
	}while(doneflag == 0 && iteration < max_iterations);

	if(doneflag == 1)
	{
		printf("Solution converged after %i iterations.\n",iteration);
		
		sprintf(momentsoutname,"%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%i%s%i%s","../data/moments_sigstim",sigstim,"_sigin",sigin,"_sigout",sigout,"_vin",vin,"_vout",vout,"_kappa",Poissk,"_veff",veff,"_seed",seed,"_iter",iteration,".txt");
		momentsoutput = fopen(momentsoutname,"w");
		fprintf(momentsoutput,"%g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n",D[1],D[2],offsets[0],offsets[1],fe_moments->favg,fi_moments->favg,fe_moments->f2avg,fi_moments->f2avg,fe_moments->sfavg, fi_moments->sfavg, fe_moments->fefiavg,fe_moments->chi2); //print moments & D's.
		fclose(momentsoutput);		
	}

	return 0;
}
