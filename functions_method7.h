// Set of functions to be used in solution of Optimal Linear Estimator solution of neuronal nonlinearity.
// Includes probability densities and integrals.

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>

struct rootfindingparams
{
	gsl_spline *splineptr; // Pointer to the interpolating spline.
	gsl_interp_accel *accptr; // Pointer to the interpolant accelerator table.
	double root; // y-value for which y = interpolant(x).
};

struct pdf_params
{
	double stimavg; // average stimulus (usually 0)
	double sigstim; // stimulus variance
	double avgin; // average input noise (usually 0)
	double sigin; // input noise variance
	double vin; // input noise correlation coefficient
	double avgout; // average output noise (usually 0)
	double sigout; // output noise variance
	double vout; // output noise correlation coefficient.
	double Poissk; // Poisson response size
	double z; // Variable argument, to be used in integration kernels.
};

struct interpolantintegrationparams
{
	gsl_integration_glfixed_table * GLtable; // Pointer to the Gauss-Legendre table.
	gsl_spline *splineptr; // Pointer to the interpolating spline.
	gsl_interp_accel *accptr; // Pointer to the interpolant accelerator table.
	int numintpts; // numpts*xcount gives number of points to use in integration table
	double zeroroot; // Value at which interpolating function crosses 0.
	double unityroot; // Value at which interpolating function crosses 1.
	int xcount; // Number of points in xdata
	struct pdf_params modelparams;
};

struct fefi_doubleintegralparams
{
	gsl_integration_glfixed_table * GLtable_fe; // Pointer to the Gauss-Legendre table.
	gsl_spline *splineptr_fe; // Pointer to the interpolating spline.
	gsl_interp_accel *accptr_fe; // Pointer to the interpolant accelerator table.
	double zeroroot_fe; // Value at which interpolating function crosses 0.
	double unityroot_fe; // Value at which interpolating function crosses 1.
	int fe_count; // Number of points in xdata

	gsl_integration_glfixed_table * GLtable_fi; // Pointer to the Gauss-Legendre table.
	gsl_spline *splineptr_fi; // Pointer to the interpolating spline.
	gsl_interp_accel *accptr_fi; // Pointer to the interpolant accelerator table.
	double zeroroot_fi; // Value at which interpolating function crosses 0.
	double unityroot_fi; // Value at which interpolating function crosses 1.
	int fi_count; // Number of points in xdata
	int numintpts; // numpts*xcount gives number of points to use in integration table
	struct pdf_params modelparams;
};

double heaviside(double x)
{
	// Heaviside step function \Theta(x). We'll use the convention that
	// \Theta(0) = 0.

	if(x > 0)
		return 1.0;
	else
		return 0.0;
}

double maxval(double a, double b)
{
	// Returns maximum of a and b.

	if(a > b)
		return a;
	else
		return b;
}

double minval(double a, double b)
{
	// Returns minimum of a and b.

	if(a < b)
		return a;
	else
		return b;
}

double rangethreshold(double x, double thresh)
{

	// This function restricts an input to the range (0,1) by
	// mapping all points > 1 to 1 and all points < 0 to 0.

	if(x >= thresh)
		return thresh;
	else
	{
		if(x <= 0)
			return 0.0;
		else
			return x;
	}
	
}

double fe_constrained(double v, void *params)
{
	// Compute fe(z) with constrains 0 <= fe(z) <= 1 imposed.
	
	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fe = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fe = intparams_fe->GLtable;
	gsl_spline * fe_spline = intparams_fe->splineptr;
	gsl_interp_accel * fe_acc = intparams_fe->accptr;
	int numpts = intparams_fe->numintpts;
	double ze_0 = intparams_fe->zeroroot;
	double ze_1 = intparams_fe->unityroot;
	int fecount = intparams_fe->xcount;
	struct pdf_params modelparams = intparams_fe->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	double z = modelparams.z;

	double temp_fe = 0;
	if(ze_0 <= ze_1)
	{
		if(v < ze_1)
		{
			if(v > ze_0)
				temp_fe = rangethreshold(gsl_spline_eval(fe_spline, v, fe_acc),1.0);
			else
				temp_fe = 0;
		}
		else
			temp_fe = 1;
	}
	else
	{
		if(v < ze_0)
		{
			if(v > ze_1)
				temp_fe = rangethreshold(gsl_spline_eval(fe_spline, v, fe_acc),1.0);
			else
				temp_fe = 1;
		}
		else
			temp_fe = 0;
	
	}

	return temp_fe;
}

double fi_constrained(double v, void *params)
{
	// Compute fi(z) with constrains 0 <= fi(z) <= 1 imposed.
	
	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fi = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fi = intparams_fi->GLtable;
	gsl_spline * fi_spline = intparams_fi->splineptr;
	gsl_interp_accel * fi_acc = intparams_fi->accptr;
	int numpts = intparams_fi->numintpts;
	double zi_0 = intparams_fi->zeroroot;
	double zi_1 = intparams_fi->unityroot;
	int ficount = intparams_fi->xcount;
	struct pdf_params modelparams = intparams_fi->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	double z = modelparams.z;

	double temp_fi = 0;
	if(zi_0 <= zi_1)
	{	
		if(v < zi_1)
		{
			if(v > zi_0)
				temp_fi = rangethreshold(gsl_spline_eval(fi_spline, v, fi_acc),1.0);
			else
				temp_fi = 0;
		}
		else
			temp_fi = 1;
	}
	else
	{
		if(v < zi_0)
		{
			if(v > zi_1)
				temp_fi = rangethreshold(gsl_spline_eval(fi_spline, v, fi_acc),1.0);
			else
				temp_fi = 1;
		}
		else
			temp_fi = 0;	
	}

	return temp_fi;
}

double stimilusdensity_gaussian(double v, void *params)
{

	// Define 
	double *p = (double *) params;
	double x = v;

	// Label parameters

	double avg = p[0];
	double var = p[1]*p[1];

	double q = exp(-(x-avg)*(x-avg)/2/var)/sqrt(2*PI*var);

	return q; 
}

double prenoisedensity_gaussian(double * v, void *params)
{

	// Define 
	double *p = (double *) params;
	double *x = v;

	// Label parameters

	double avg0 = p[0];
	double var0 = p[1]*p[1];
	double avg1 = p[2];
	double var1 = p[3]*p[3];
	double rho = p[4]; // correlation coefficient
	double hyprho = 1 - rho*rho; // convenient grouping

	double q = exp(-((x[0]-avg0)*(x[0]-avg0)/var0 + (x[1]-avg1)*(x[1]-avg1)/var1 - 2*rho*(x[0]-avg0)*(x[1]-avg1)/sqrt(var0*var1))/2/hyprho)/(2*PI)/sqrt(var0*var1*hyprho);

	return q; 
}

double postnoisedensity_gaussian(double * v, void *params)
{

	// Define 
	double *p = (double *) params;
	double *x = v;

	// Label parameters

	double avg0 = p[0];
	double var0 = p[1]*p[1];
	double avg1 = p[2];
	double var1 = p[3]*p[3];
	double rho = p[4]; // correlation coefficient
	double hyprho = 1 - rho*rho; // convenient grouping

	double q = exp(-((x[0]-avg0)*(x[0]-avg0)/var0 + (x[1]-avg1)*(x[1]-avg1)/var1 - 2*rho*(x[0]-avg0)*(x[1]-avg1)/sqrt(var0*var1))/2/hyprho)/(2*PI)/sqrt(var0*var1*hyprho);

	return q; 
}

double g_gaussian(double *v, void * params)
{
	// function g

	double x = v[0];
	double y = v[1];

	double *p = (double *) params;
	double sigstim = p[0], sigin = p[1], vin = p[2];

	double argnum = sigstim*sigstim*(x+y)*(x+y) + sigin*sigin*(x*x+y*y-2.0*vin*x*y);
	double argden = (1.0+vin)*sigin*sigin*(2.0*sigstim*sigstim+(1.0-vin)*sigin*sigin);

	double g = exp(-0.5*argnum/argden)/2.0/PI/sqrt(argden);

	return g;
}

double K1_gaussiankernel(double v, void * params)
{
	// Integral kernel K1

	// Doesn't reuse function for g as in needs to treat x as a parameter
	// (this function K1 is for use with the gsl integrator).
	
		// **THIS FUNCTION HAS BEEN CORRECTED WITH THE CORRECT SIGNS

	double y = v;

	struct pdf_params * modelparams = (struct pdf_params *) params;

	double stimavg=modelparams->stimavg, sigstim = modelparams->sigstim; 
	double avgin=modelparams->avgin, sigin = modelparams->sigin, vin = modelparams->vin; 
	double avgout=modelparams->avgout, sigout = modelparams->sigout, vout = modelparams->vout;
	double x = modelparams->z;
	double svar = sigstim*sigstim, invar = sigin*sigin;

	double argnum = (svar+invar)*(x*x+y*y) -2*(svar+invar*vin)*x*y;
	double argden = (1.0-vin)*invar*(2.0*svar+(1.0+vin)*invar);

	double g = exp(-0.5*argnum/argden)/2.0/PI/sqrt(argden);

	double K = sqrt(2*PI)*sqrt(svar+invar)*exp(0.5*x*x/(svar+invar))*g;

	return K;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Analytically computed integrals over the kernel K1 with closed forms
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double upperintegratedK1_gaussiankernel(double v, void * params)
{
	// Integral of the kernel K1(x,z) from -\infty to a.
	// = Integral of kernel K1(-x,z) from -a to \infty
	// The integral has an exact representation as a special function,
	// so we use that here to limit integration error.
	
	// **THIS FUNCTION HAS BEEN CORRECTED WITH THE CORRECT SIGNS


	double a = -v; // since the function is passed z1, but a is -z1.

	struct pdf_params * modelparams = (struct pdf_params *) params;

	double stimavg=modelparams->stimavg, sigstim = modelparams->sigstim; 
	double avgin=modelparams->avgin, sigin = modelparams->sigin, vin = modelparams->vin; 
	double avgout=modelparams->avgout, sigout = modelparams->sigout, vout = modelparams->vout;
	double z = modelparams->z;
	double svar = sigstim*sigstim, invar = sigin*sigin;

	double argnum = a*(svar+invar)-(svar+invar*vin)*z;
	double argden = sqrt((2.0*svar+(1.0+vin)*invar)*(1.0-vin)*2.0*invar*(svar+invar));

	double err = gsl_sf_erf(argnum/argden);

	double integral = 0.5*(1+err);

	return integral;
}

// The following function is not presently necessary

//double lowerintegratedK1_gaussiankernel(double v, void * params)
//{
	// Integral of the kernel K1(x,z) from a to +\infty.
	// The integral has an exact representation as a special function,
	// so we use that here to limit integration error.


//	double a = v;

//	struct pdf_params * modelparams = (struct pdf_params *) params;

//	double stimavg=modelparams->stimavg, sigstim = modelparams->sigstim; 
//	double avgin=modelparams->avgin, sigin = modelparams->sigin, vin = modelparams->vin; 
//	double avgout=modelparams->avgout, sigout = modelparams->sigout, vout = modelparams->vout;
//	double z = modelparams->z;

//	double argnum = a*(sigstim*sigstim+sigin*sigin)+(sigstim*sigstim-sigin*sigin*vin)*z;
//	double argden = sqrt((2.0*sigstim*sigstim+(1.0-vin)*sigin*sigin)*(1.0+vin))*2.0*sigin;

//	double err = gsl_sf_erf(argnum/argden);

//	double integral = 0.5*(1-err);

//	return integral;
//}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Functions involving interpolants. Most functions are integrals over the interpolants
//	and the kernel K1(x,y) or the function g(x,y) using Gauss-Legendre integration.
//
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

double interpolantminusxval(double x, void * params)
{

	// This function is to be used in a root-finding solver. This function
	// is passed an argument x, which is the value at which it is to be evaluated,
	// and a struct containing pointers to the gsl interpolant functions, as well
	// as the desired y value for which we want to solve y = interpolant(x).

	// For convenience, this is the definition of the struct:

	// struct rootfindingparams
	// {
	//	gsl_spline *splineptr; // GSL spline pointer
	//	gsl_interp_accel *accptr; // GSL spline accelerator pointer
	//	double root; // Desired y-value, y = interpolant(x).
	//};

	// See http://www.gnu.org/software/gsl/manual/html_node/Interpolation.html for more info.

	struct rootfindingparams *p = (struct rootfindingparams *) params;

	double yval = gsl_spline_eval(p->splineptr, x, p->accptr)-p->root;

	return yval;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fe_avg(void * params)
{
	// This function calculates the average <f_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < f_E(s+\eta_E) > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fe(x)

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fe = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fe = intparams_fe->GLtable;
	gsl_spline * fe_spline = intparams_fe->splineptr;
	gsl_interp_accel * fe_acc = intparams_fe->accptr;
	int numpts = intparams_fe->numintpts;
	double ze_0 = intparams_fe->zeroroot;
	double ze_1 = intparams_fe->unityroot;
	int fecount = intparams_fe->xcount;
	struct pdf_params modelparams = intparams_fe->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	
	// Define some parameter combinations
	double scale = sqrt(2.0*(sigstim*sigstim + sigin*sigin));

	// Define some internal variables.

	int fe_i=0; // loop counter
	double fe_sum = 0; // Sum to approximate integral.
	double fe_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double temp_fe = 0; // Hold value of fe(z)
	
	double zlow = -4; //minval(ze_0,ze_1);
	double zhigh = 4; //maxval(ze_0,ze_1);


	for(fe_i=0;fe_i<numpts*fecount;fe_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, fe_i, &xdummy, &weightdummy, t_fe);
		// Above command returns integration at (fe2_i)th integration point xdummy, for an integral between
		// the limits of ze_0 and ze_1.

		// Integrate (sum): weight_i*interpolant(x_i)^2*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		temp_fe = fe_constrained(scale*xdummy,intparams_fe); 
		fe_sum += weightdummy*temp_fe*exp(-xdummy*xdummy);  
	}
	fe_sum = fe_sum/sqrt(PI); // multiply by constant factors.	
	
	return fe_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fi_avg(void * params)
{
	// This function calculates the average <f_I(-s+\eta_I)> using an
	// interpolant that approximates f_I(z) on its domain for which it is between
	// 0 < f_I < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_I(z).

	// The integral is calculated as: 
	// < f_I(-s+\eta_I)^2 > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * f_I(x)

	// Because 0 <= f_I(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fi = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fi = intparams_fi->GLtable;
	gsl_spline * fi_spline = intparams_fi->splineptr;
	gsl_interp_accel * fi_acc = intparams_fi->accptr;
	int numpts = intparams_fi->numintpts;
	double zi_0 = intparams_fi->zeroroot;
	double zi_1 = intparams_fi->unityroot;
	int ficount = intparams_fi->xcount;
	struct pdf_params modelparams = intparams_fi->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	
	// Define some parameter combinations
	double scale = sqrt(2.0*(sigstim*sigstim + sigin*sigin));

	// Define some internal variables.

	int fi_i=0; // loop counter
	double fi_sum = 0; // Sum to approximate integral.
	double fi_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_I(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	double temp_fi = 0; // Hold fi(z) value.
	
	double zlow = -4; //minval(zi_0,zi_1);
	double zhigh = 4; //maxval(zi_0,zi_1);


	for(fi_i=0;fi_i<numpts*ficount;fi_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, fi_i, &xdummy, &weightdummy, t_fi);
		// Above command returns integration at (sfe_i)th integration point xdummy, for an integral between
		// the limits of ze_data[0] and ze_data[fecount-1].

		// Integrate (sum): weight_i*interpolant(x_i)*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))

		temp_fi = fi_constrained(scale*xdummy,intparams_fi);
		fi_sum += weightdummy*temp_fi*exp(-xdummy*xdummy); 
	}
	fi_sum = fi_sum/sqrt(PI); // Multiply by constant factors.

	return fi_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double sfe_avg(void * params)
{
	// This function calculates the average <s f_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < s f_E(s+\eta_E) > = \sigma_s^2/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{3/2} \times
	// \int_{-\infty}^\infty dx x \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * f_E(x)

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fe = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fe = intparams_fe->GLtable;
	gsl_spline * fe_spline = intparams_fe->splineptr;
	gsl_interp_accel * fe_acc = intparams_fe->accptr;
	int numpts = intparams_fe->numintpts;
	double ze_0 = intparams_fe->zeroroot;
	double ze_1 = intparams_fe->unityroot;
	int fecount = intparams_fe->xcount;
	struct pdf_params modelparams = intparams_fe->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	
	// Define some parameter combinations
	double scale = sqrt(2.0*(sigstim*sigstim + sigin*sigin));

	// Define some internal variables.

	int sfe_i=0; // loop counter
	double sfe_sum = 0; // Sum to approximate integral.
	double sfe_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double temp_fe = 0; // Hold fe(z)
	
	double zlow = -4; //minval(ze_0,ze_1);
	double zhigh = 4; //maxval(ze_0,ze_1);

	for(sfe_i=0;sfe_i<numpts*fecount;sfe_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, sfe_i, &xdummy, &weightdummy, t_fe);
		// Above command returns integration at (sfe_i)th integration point xdummy, for an integral between
		// the limits of ze_data[0] and ze_data[fecount-1].

		// Integrate (sum): weight_i*interpolant(x_i)*x_i*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		temp_fe = fe_constrained(scale*xdummy,intparams_fe);
		sfe_sum += weightdummy*temp_fe*xdummy*exp(-xdummy*xdummy);  
	}
	
	sfe_sum = 2.0*sigstim*sigstim/scale/sqrt(PI)*sfe_sum; // Need to multiply by factor out front of integral.	

	return sfe_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double sfi_avg(void * params)
{
	// This function calculates the average <s f_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < s f_E(s+\eta_E) > = \sigma_s^2/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{3/2} \times
	// \int_{-\infty}^\infty dx x \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * f_E(x)

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fi = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fi = intparams_fi->GLtable;
	gsl_spline * fi_spline = intparams_fi->splineptr;
	gsl_interp_accel * fi_acc = intparams_fi->accptr;
	int numpts = intparams_fi->numintpts;
	double zi_0 = intparams_fi->zeroroot;
	double zi_1 = intparams_fi->unityroot;
	int ficount = intparams_fi->xcount;
	struct pdf_params modelparams = intparams_fi->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

	// Define some internal variables.

	int sfi_i=0; // loop counter
	double sfi_sum = 0; // Sum to approximate integral.
	double sfi_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double temp_fi = 0; // Hold fi(z)
	
	// Define some parameter combinations
	double scale = sqrt(2.0*(sigstim*sigstim + sigin*sigin));
	
	double zlow = -4; //minval(zi_0,zi_1);
	double zhigh = 4; //maxval(zi_0,zi_1);

	for(sfi_i=0;sfi_i<numpts*ficount;sfi_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, sfi_i, &xdummy, &weightdummy, t_fi);
		// Above command returns integration at (sfe_i)th integration point xdummy, for an integral between
		// the limits of ze_data[0] and ze_data[fecount-1].

		// Integrate (sum): weight_i*interpolant(x_i)*x_i*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		temp_fi = fi_constrained(scale*xdummy,intparams_fi);
		sfi_sum += weightdummy*temp_fi*xdummy*exp(-xdummy*xdummy);  
	}
	
	sfi_sum = -2.0*sigstim*sigstim/scale/sqrt(PI)*sfi_sum; // Need to multiply by factor out front of integral.	

	return sfi_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fe2_avg(void * params)
{
	// This function calculates the average <f^2_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < f^2_E(s+\eta_E) > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fe(x)^2

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fe = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fe = intparams_fe->GLtable;
	gsl_spline * fe_spline = intparams_fe->splineptr;
	gsl_interp_accel * fe_acc = intparams_fe->accptr;
	int numpts = intparams_fe->numintpts;
	double ze_0 = intparams_fe->zeroroot;
	double ze_1 = intparams_fe->unityroot;
	int fecount = intparams_fe->xcount;
	struct pdf_params modelparams = intparams_fe->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	
	// Define some parameter combinations
	double scale = sqrt(2.0*(sigstim*sigstim + sigin*sigin));

	// Define some internal variables.

	int fe2_i=0; // loop counter
	double fe2_sum = 0; // Sum to approximate integral.
	double fe2_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	double temp_fe = 0; // Hold fe(z)
	
	
	double zlow = -4; //minval(ze_0,ze_1);
	double zhigh = 4; //maxval(ze_0,ze_1);


	for(fe2_i=0;fe2_i<numpts*fecount;fe2_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, fe2_i, &xdummy, &weightdummy, t_fe);
		// Above command returns integration at (fe2_i)th integration point xdummy, for an integral between
		// the limits of ze_0 and ze_1.

		temp_fe = fe_constrained(scale*xdummy,intparams_fe);
		fe2_sum += weightdummy*temp_fe*temp_fe*exp(-xdummy*xdummy);  
	}
	fe2_sum = fe2_sum/sqrt(PI); // multiply by constant factors.	

	return fe2_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fi2_avg(void * params)
{
	// This function calculates the average <f^2_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < f^2_E(s+\eta_E) > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fe(x)^2

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fi = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fi = intparams_fi->GLtable;
	gsl_spline * fi_spline = intparams_fi->splineptr;
	gsl_interp_accel * fi_acc = intparams_fi->accptr;
	int numpts = intparams_fi->numintpts;
	double zi_0 = intparams_fi->zeroroot;
	double zi_1 = intparams_fi->unityroot;
	int ficount = intparams_fi->xcount;
	struct pdf_params modelparams = intparams_fi->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	
	// Define some parameter combinations
	double scale = sqrt(2.0*(sigstim*sigstim + sigin*sigin));

	// Define some internal variables.

	int fi2_i=0; // loop counter
	double fi2_sum = 0; // Sum to approximate integral.
	double fi2_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	double temp_fi = 0; // Hold fi(z)
	
	double zlow = -4; //minval(zi_0,zi_1);
	double zhigh = 4; //maxval(zi_0,zi_1);

	for(fi2_i=0;fi2_i<numpts*ficount;fi2_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, fi2_i, &xdummy, &weightdummy, t_fi);
		// Above command returns integration at (fe2_i)th integration point xdummy, for an integral between
		// the limits of zi_0 and zi_1.

		// Integrate (sum): weight_i*interpolant(x_i)^2*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		temp_fi = fi_constrained(scale*xdummy,intparams_fi);
		fi2_sum += weightdummy*temp_fi*temp_fi*exp(-xdummy*xdummy);  
	}
	fi2_sum = fi2_sum/sqrt(PI); // multiply by constant factors.	
	
	return fi2_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fegfi_montecarlo_integrand(double *k, size_t dim, void *params)
{

	// This is the integral of the two-dimensional integral we need to perform
	// to calculate the average <f_E(s+\eta_E) f_I(-s-\eta_I)>. We have changed variables
	// to x = (\sqrt{2}*u - |b|*w)/2 and y = (\sqrt{2}*u + |b|*w)/2 in order to factorize the
	// function g(x,y) = h(u)h(w). Though the terms f_E(x) and f_I(y) are no longer functions of
	// single variables, this shouldn't present any problems for the integration.

	// Get parameters and define convenient names:
	
	struct fefi_doubleintegralparams * doubleparams = (struct fefi_doubleintegralparams *) params;

	gsl_integration_glfixed_table * t_fi = doubleparams->GLtable_fi;
	gsl_spline * fi_spline = doubleparams->splineptr_fi;
	gsl_interp_accel * fi_acc = doubleparams->accptr_fi;
	double zi_0 = doubleparams->zeroroot_fi;
	double zi_1 = doubleparams->unityroot_fi;
	int ficount = doubleparams->fi_count;

	gsl_integration_glfixed_table * t_fe = doubleparams->GLtable_fe;
	gsl_spline * fe_spline = doubleparams->splineptr_fe;
	gsl_interp_accel * fe_acc = doubleparams->accptr_fe;
	double ze_0 = doubleparams->zeroroot_fe;
	double ze_1 = doubleparams->unityroot_fe;
	int fecount = doubleparams->fe_count;
	
	int numpts = doubleparams->numintpts;
	struct pdf_params modelparams = doubleparams->modelparams;

	//	printf("monte carlo integrand func check: z0 = %g, z1 = %g\n",ze_0,ze_1);

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

	double scale = sqrt(2.0*(sigstim*sigstim+sigin*sigin));
	double nueff = (sigstim*sigstim+sigin*sigin*vin)/(sigstim*sigstim+sigin*sigin);

	double argument_plus = scale*(sqrt((1.0-nueff)/2.0)*k[0]+sqrt((1.0+nueff)/2.0)*k[1]);
	double argument_minus = scale*(sqrt((1.0-nueff)/2.0)*k[0]-sqrt((1.0+nueff)/2.0)*k[1]);
	
//	printf("arg+ = %g, arg- = %g\n",argument_plus,argument_minus);
	
	double temp_fe = 0; 
	double temp_fi = 0;
	
	// fe_constrained and fi_constrained take interpolantintegrationparams struct, not 
	// doubleparams struct, so define appropriate structs to pass.
	struct interpolantintegrationparams params_fe;
	params_fe.GLtable = t_fe;
	params_fe.splineptr = fe_spline;
	params_fe.accptr = fe_acc;
	params_fe.numintpts = numpts;
	params_fe.zeroroot = ze_0;
	params_fe.unityroot = ze_1;
	params_fe.xcount = fecount;
	params_fe.modelparams = modelparams;
	
	struct interpolantintegrationparams params_fi;
	params_fi.GLtable = t_fi;
	params_fi.splineptr = fi_spline;
	params_fi.accptr = fi_acc;
	params_fi.numintpts = numpts;
	params_fi.zeroroot = zi_0;
	params_fi.unityroot = zi_1;
	params_fi.xcount = ficount;
	params_fi.modelparams = modelparams;
	
	temp_fe = fe_constrained(argument_plus,&params_fe);
	temp_fi = fi_constrained(argument_minus,&params_fi);
			
	double A = temp_fe*temp_fi*exp(-(k[0]*k[0]+k[1]*k[1]))/PI; 

	return A;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fegfi_montecarlo_integral(void * params)
{

	// This function calculates the average <f_E(s+\eta_E) f_I(-s-\eta_I)> using an
	// interpolant that approximates f_E(z) and f_I(z) on the domains for which they are between
	// 0 < f < 1. The integration uses a plain Monte Carlo integration scheme (rather than an adaptive one
	// So that the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E/I(z).

	// The integral is a double integral, equal to : 
	// < f_E(s+\eta_E) f_I(-s+\eta_I) > = 
	// \int_{-\infty}^\infty dx dy f_E(x)*g(x,y) * f_I(y)
	// where g(x,y) is a function defined earlier. This is a double integral, hence the use of a Monte Carlo
	// procedure for calculating it. Because g(x,y) is a multi-Gaussian, it will fall off to nearly zero at some x or y. We can
	// cut off our integration regions at these points. However, for different x (for example), the kernel does
	// not always decay to ~0 at the same values of y. To be careful about it, we will actually perform a
	// change of variables (x,y) -> (u,w) which factorizes the integral kernel into separate functions of
	// u and w, both the same Gaussian. In terms of the u and w infinite limits, we can easily cut off the
	// integrals at fixed values.

	// Under the change of variables \sqrt{2} u = x + y, w = (y - x)/|b|, 
	// where |b| = \sqrt{2*(2*sigma_s^2+sigma_{in}^2*(1+\nu_{in}))/\sigma_{in}^2/(1-\nu_{in})},
	// we have dx dy g(x,y) -> dw du h(u)h(w), where we have absorbed the Jacobian factor into the definitions of
	// h(.), which has the Gaussian form
	// h(x) = \exp(-1/2 x^2/(\sigma_{in}^2(1-\nu_{in})). 99% of the area under this curve is contained with 
	// x = \pm 3 \sigma_{in} \sqrt{1-\nu_{in}), so we will replace all infinite limits over u or w with this limit.

	// Get parameters and define convenient names:

	struct fefi_doubleintegralparams * doubleparams = (struct fefi_doubleintegralparams *) params;

	gsl_integration_glfixed_table * t_fi = doubleparams->GLtable_fi;
	gsl_spline * fi_spline = doubleparams->splineptr_fi;
	gsl_interp_accel * fi_acc = doubleparams->accptr_fi;
	double zi_0 = doubleparams->zeroroot_fi;
	double zi_1 = doubleparams->unityroot_fi;
	int ficount = doubleparams->fi_count;

	gsl_integration_glfixed_table * t_fe = doubleparams->GLtable_fe;
	gsl_spline * fe_spline = doubleparams->splineptr_fe;
	gsl_interp_accel * fe_acc = doubleparams->accptr_fe;
	double ze_0 = doubleparams->zeroroot_fe;
	double ze_1 = doubleparams->unityroot_fe;
	int fecount = doubleparams->fe_count;
	
	int numpts = doubleparams->numintpts;
	struct pdf_params modelparams = doubleparams->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

  	double res, err;

  	double xl[2] = { -4, -4};
  	double xu[2] = { 4, 4};

  	const gsl_rng_type *T;
  	gsl_rng *r;

  	gsl_monte_function G = { &fegfi_montecarlo_integrand, 2, doubleparams };

  	size_t calls = 10*1000000;

  	gsl_rng_env_setup ();

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	  
	gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
	gsl_monte_plain_integrate (&G, xl, xu, 2, calls, r, s,&res, &err);
	gsl_monte_plain_free (s);

	gsl_rng_free (r);

	//printf("monte carlo integral result = %g, error = %g\n",res,err);

	return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double feK1integral(void * params)
{
	// This function calculates the average <s f_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < s f_E(s+\eta_E) > = \sigma_s^2/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{3/2} \times
	// \int_{-\infty}^\infty dx x \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * f_E(x)

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fe = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fe = intparams_fe->GLtable;
	gsl_spline * fe_spline = intparams_fe->splineptr;
	gsl_interp_accel * fe_acc = intparams_fe->accptr;
	int numpts = intparams_fe->numintpts;
	double ze_0 = intparams_fe->zeroroot;
	double ze_1 = intparams_fe->unityroot;
	int fecount = intparams_fe->xcount;
	struct pdf_params modelparams = intparams_fe->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	double z = modelparams.z;
	
	// Parameter combinations
	double scale = sqrt(2.0*(sigstim*sigstim + sigin*sigin));
	double nueff = (sigstim*sigstim + sigin*sigin*vin)/(sigstim*sigstim + sigin*sigin);

	// Define some internal variables.

	int feK1_i=0; // loop counter
	double feK1_sum = 0; // Sum to approximate integral.
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.

	double zlow = -4; //minval(ze_0,ze_1);
	double zhigh = 4; //maxval(ze_0,ze_1);
	
	double temp_fe=0; // hold value of f_e(z)
	double argument=0; // hold value of argument to fe(z)

	for(feK1_i=0;feK1_i<numpts*fecount;feK1_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, feK1_i, &xdummy, &weightdummy, t_fe);
		// Above command returns integration at (feK1_i)th integration point xdummy, for an integral between
		// the limits of zlow and zhigh.
		
		argument = scale*sqrt(1.0-nueff*nueff)*xdummy - nueff*z;
		temp_fe = fe_constrained(argument,intparams_fe);

		// Integrate (sum): weight_i*interpolant(x_i)*K1(-x_i,z)
		feK1_sum += weightdummy*temp_fe*exp(-xdummy*xdummy);  
	}
	
	feK1_sum = feK1_sum/sqrt(PI);

	return feK1_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fiK1integral(void * params)
{
	// This function calculates the average <s f_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < s f_E(s+\eta_E) > = \sigma_s^2/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{3/2} \times
	// \int_{-\infty}^\infty dx x \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * f_E(x)

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fi = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fi = intparams_fi->GLtable;
	gsl_spline * fi_spline = intparams_fi->splineptr;
	gsl_interp_accel * fi_acc = intparams_fi->accptr;
	int numpts = intparams_fi->numintpts;
	double zi_0 = intparams_fi->zeroroot;
	double zi_1 = intparams_fi->unityroot;
	int ficount = intparams_fi->xcount;
	struct pdf_params modelparams = intparams_fi->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;
	double z = modelparams.z;
	
	// Parameter combinations
	double scale = sqrt(2.0*(sigstim*sigstim + sigin*sigin));
	double nueff = (sigstim*sigstim + sigin*sigin*vin)/(sigstim*sigstim + sigin*sigin);

	// Define some internal variables.

	int fiK1_i=0; // loop counter
	double fiK1_sum = 0; // Sum to approximate integral.
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double zlow = -4; //minval(zi_0,zi_1);
	double zhigh = 4; //maxval(zi_0,zi_1);

	double temp_fi=0; // hold value of fi(z)
	double argument=0; // hold value of argument to fi(z)
	for(fiK1_i=0;fiK1_i<numpts*ficount;fiK1_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, fiK1_i, &xdummy, &weightdummy, t_fi);
		// Above command returns integration at (feK1_i)th integration point xdummy, for an integral between
		// the limits of ze_0 and ze_1.
		
		argument = scale*sqrt(1.0-nueff*nueff)*xdummy - nueff*z;
		temp_fi = fi_constrained(argument,intparams_fi);

		// Integrate (sum): weight_i*interpolant(x_i)*K1(-x_i,z)
		fiK1_sum += weightdummy*temp_fi*exp(-xdummy*xdummy);  
	}
	
	fiK1_sum = fiK1_sum/sqrt(PI);

	return fiK1_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
double fe_avg_old(void * params)
{
	// This function calculates the average <f_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < f_E(s+\eta_E) > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fe(x)

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fe = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fe = intparams_fe->GLtable;
	gsl_spline * fe_spline = intparams_fe->splineptr;
	gsl_interp_accel * fe_acc = intparams_fe->accptr;
	int numpts = intparams_fe->numintpts;
	double ze_0 = intparams_fe->zeroroot;
	double ze_1 = intparams_fe->unityroot;
	int fecount = intparams_fe->xcount;
	struct pdf_params modelparams = intparams_fe->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

	// Define some internal variables.

	int fe_i=0; // loop counter
	double fe_sum = 0; // Sum to approximate integral.
	double fe_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double zlow = minval(ze_0,ze_1);
	double zhigh = maxval(ze_0,ze_1);

	// Computing the term \int_{z_E1}^\infty dx exp[-x^2/2(\sigma_s^2+\sigma_{in}^2)]/sqrt(2*pi*(sigstim^2+sigin^2)) (closed form expression).
	fe_constpart = 0.5*(1-gsl_sf_erf(ze_1/sqrt(2*(sigstim*sigstim+sigin*sigin)))); // (1/2)*(1 - erf(a/sqrt(2*(sigstim^2+sigin^2))))

	for(fe_i=0;fe_i<numpts*fecount;fe_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, fe_i, &xdummy, &weightdummy, t_fe);
		// Above command returns integration at (fe2_i)th integration point xdummy, for an integral between
		// the limits of ze_0 and ze_1.

		// Integrate (sum): weight_i*interpolant(x_i)^2*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		double temp_fe = gsl_spline_eval(fe_spline, xdummy, fe_acc);
		fe_sum += weightdummy*rangethreshold(temp_fe,1.0)*exp(-0.5*xdummy*xdummy/(sigstim*sigstim+sigin*sigin));  
	}
	fe_sum = fe_sum/sqrt(2*PI*(sigstim*sigstim+sigin*sigin)); // multiply by constant factors.	
	
	if(ze_0 <= ze_1)
		fe_sum += fe_constpart;
	else
		fe_sum += 1-fe_constpart;

	return fe_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fi_avg_old(void * params)
{
	// This function calculates the average <f_I(-s+\eta_I)> using an
	// interpolant that approximates f_I(z) on its domain for which it is between
	// 0 < f_I < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_I(z).

	// The integral is calculated as: 
	// < f_I(-s+\eta_I)^2 > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * f_I(x)

	// Because 0 <= f_I(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fi = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fi = intparams_fi->GLtable;
	gsl_spline * fi_spline = intparams_fi->splineptr;
	gsl_interp_accel * fi_acc = intparams_fi->accptr;
	int numpts = intparams_fi->numintpts;
	double zi_0 = intparams_fi->zeroroot;
	double zi_1 = intparams_fi->unityroot;
	int ficount = intparams_fi->xcount;
	struct pdf_params modelparams = intparams_fi->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

	// Define some internal variables.

	int fi_i=0; // loop counter
	double fi_sum = 0; // Sum to approximate integral.
	double fi_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_I(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double zlow = minval(zi_0,zi_1);
	double zhigh = maxval(zi_0,zi_1);

	// Computing the term \int_{-\infty}^{zI1} dx exp[-x^2/2(\sigma_s^2+\sigma_{in}^2)]/sqrt(2*pi*(sigstim^2+sigin^2)) (closed form expression).

	fi_constpart = 0.5*(1-gsl_sf_erf(zi_1/sqrt(2*(sigstim*sigstim+sigin*sigin)))); // (1/2)*(1 + erf(a/sqrt(2*(sigstim^2+sigin^2))))

	for(fi_i=0;fi_i<numpts*ficount;fi_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, fi_i, &xdummy, &weightdummy, t_fi);
		// Above command returns integration at (sfe_i)th integration point xdummy, for an integral between
		// the limits of ze_data[0] and ze_data[fecount-1].

		// Integrate (sum): weight_i*interpolant(x_i)*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		double temp_fi = gsl_spline_eval(fi_spline, xdummy, fi_acc);
		fi_sum += weightdummy*rangethreshold(temp_fi,1.0)*exp(-0.5*xdummy*xdummy/(sigstim*sigstim+sigin*sigin)); 
	}
	fi_sum = fi_sum/sqrt(2*PI*(sigstim*sigstim+sigin*sigin)); // Multiply by constant factors.

	if(zi_0 <= zi_1)
		fi_sum += fi_constpart;
	else
		fi_sum += 1-fi_constpart;

	return fi_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double sfe_avg_old(void * params)
{
	// This function calculates the average <s f_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < s f_E(s+\eta_E) > = \sigma_s^2/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{3/2} \times
	// \int_{-\infty}^\infty dx x \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * f_E(x)

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fe = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fe = intparams_fe->GLtable;
	gsl_spline * fe_spline = intparams_fe->splineptr;
	gsl_interp_accel * fe_acc = intparams_fe->accptr;
	int numpts = intparams_fe->numintpts;
	double ze_0 = intparams_fe->zeroroot;
	double ze_1 = intparams_fe->unityroot;
	int fecount = intparams_fe->xcount;
	struct pdf_params modelparams = intparams_fe->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

	// Define some internal variables.

	int sfe_i=0; // loop counter
	double sfe_sum = 0; // Sum to approximate integral.
	double sfe_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double zlow = minval(ze_0,ze_1);
	double zhigh = maxval(ze_0,ze_1);

	// Computing the term \int_{z_E1}^\infty dx x exp[-x^2/2(\sigma_s^2+\sigma_{in}^2)] (closed form expression).
	sfe_constpart = sigstim*sigstim/sqrt(2*PI*(sigstim*sigstim+sigin*sigin))*exp(-ze_1*ze_1/2/(sigstim*sigstim+sigin*sigin));

	for(sfe_i=0;sfe_i<numpts*fecount;sfe_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, sfe_i, &xdummy, &weightdummy, t_fe);
		// Above command returns integration at (sfe_i)th integration point xdummy, for an integral between
		// the limits of ze_data[0] and ze_data[fecount-1].

		// Integrate (sum): weight_i*interpolant(x_i)*x_i*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		double temp_fe = gsl_spline_eval(fe_spline, xdummy, fe_acc);
		sfe_sum += weightdummy*rangethreshold(temp_fe,1.0)*xdummy*exp(-0.5*xdummy*xdummy/(sigstim*sigstim+sigin*sigin));  
	}
	
	sfe_sum = sigstim*sigstim/sqrt(2*PI*(sigstim*sigstim+sigin*sigin))/(sigstim*sigstim+sigin*sigin)*sfe_sum; // Need to multiply by factor out front of integral.	

	if(ze_0 <= ze_1)
		sfe_sum += sfe_constpart;
	else
		sfe_sum += -sfe_constpart;

	return sfe_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double sfi_avg_old(void * params)
{
	// This function calculates the average <s f_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < s f_E(s+\eta_E) > = \sigma_s^2/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{3/2} \times
	// \int_{-\infty}^\infty dx x \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * f_E(x)

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fi = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fi = intparams_fi->GLtable;
	gsl_spline * fi_spline = intparams_fi->splineptr;
	gsl_interp_accel * fi_acc = intparams_fi->accptr;
	int numpts = intparams_fi->numintpts;
	double zi_0 = intparams_fi->zeroroot;
	double zi_1 = intparams_fi->unityroot;
	int ficount = intparams_fi->xcount;
	struct pdf_params modelparams = intparams_fi->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

	// Define some internal variables.

	int sfi_i=0; // loop counter
	double sfi_sum = 0; // Sum to approximate integral.
	double sfi_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double zlow = minval(zi_0,zi_1);
	double zhigh = maxval(zi_0,zi_1);

	// Computing the term \int_{z_E1}^\infty dx x exp[-x^2/2(\sigma_s^2+\sigma_{in}^2)] (closed form expression).
	sfi_constpart = sigstim*sigstim/sqrt(2*PI*(sigstim*sigstim+sigin*sigin))*exp(-zi_1*zi_1/2/(sigstim*sigstim+sigin*sigin));

	for(sfi_i=0;sfi_i<numpts*ficount;sfi_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, sfi_i, &xdummy, &weightdummy, t_fi);
		// Above command returns integration at (sfe_i)th integration point xdummy, for an integral between
		// the limits of ze_data[0] and ze_data[fecount-1].

		// Integrate (sum): weight_i*interpolant(x_i)*x_i*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		double temp_fi = gsl_spline_eval(fi_spline, xdummy, fi_acc);
		sfi_sum += weightdummy*rangethreshold(temp_fi,1.0)*xdummy*exp(-0.5*xdummy*xdummy/(sigstim*sigstim+sigin*sigin));  
	}
	
	sfi_sum = sigstim*sigstim/sqrt(2*PI*(sigstim*sigstim+sigin*sigin))/(sigstim*sigstim+sigin*sigin)*sfi_sum; // Need to multiply by factor out front of integral.	

	if(zi_0 <= zi_1)
		sfi_sum += sfi_constpart;
	else
		sfi_sum += -sfi_constpart;
	
	sfi_sum = -sfi_sum; // Need to multiply by overall minus sign

	return sfi_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fe2_avg_old(void * params)
{
	// This function calculates the average <f^2_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < f^2_E(s+\eta_E) > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fe(x)^2

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fe = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fe = intparams_fe->GLtable;
	gsl_spline * fe_spline = intparams_fe->splineptr;
	gsl_interp_accel * fe_acc = intparams_fe->accptr;
	int numpts = intparams_fe->numintpts;
	double ze_0 = intparams_fe->zeroroot;
	double ze_1 = intparams_fe->unityroot;
	int fecount = intparams_fe->xcount;
	struct pdf_params modelparams = intparams_fe->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

	// Define some internal variables.

	int fe2_i=0; // loop counter
	double fe2_sum = 0; // Sum to approximate integral.
	double fe2_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double zlow = minval(ze_0,ze_1);
	double zhigh = maxval(ze_0,ze_1);

	// Computing the term \int_{z_E1}^\infty dx exp[-x^2/2(\sigma_s^2+\sigma_{in}^2)]/sqrt(2*pi*(sigstim^2+sigin^2)) (closed form expression).
	fe2_constpart = 0.5*(1-gsl_sf_erf(ze_1/sqrt(2*(sigstim*sigstim+sigin*sigin)))); // (1/2)*(1 - erf(a/sqrt(2*(sigstim^2+sigin^2))))

	for(fe2_i=0;fe2_i<numpts*fecount;fe2_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, fe2_i, &xdummy, &weightdummy, t_fe);
		// Above command returns integration at (fe2_i)th integration point xdummy, for an integral between
		// the limits of ze_0 and ze_1.

		// Integrate (sum): weight_i*interpolant(x_i)^2*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		double temp_fe = gsl_spline_eval(fe_spline, xdummy, fe_acc);
		fe2_sum += weightdummy*rangethreshold(temp_fe,1.0)*rangethreshold(temp_fe,1.0)*exp(-0.5*xdummy*xdummy/(sigstim*sigstim+sigin*sigin));  
	}
	fe2_sum = fe2_sum/sqrt(2*PI*(sigstim*sigstim+sigin*sigin)); // multiply by constant factors.	

	if(ze_0 <= ze_1)
		fe2_sum += fe2_constpart;
	else
		fe2_sum += 1-fe2_constpart;

	return fe2_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fi2_avg_old(void * params)
{
	// This function calculates the average <f^2_E(s+\eta_e)> using an
	// interpolant that approximates f_E(z) on its domain for which it is between
	// 0 < f_E < 1. The integration uses a Gauss-Legendre quadrature scheme in order
	// to match the integration points with the number of interpolant points. (So that
	// the integral is not performed at too fine a scale, at which interpolation features
	// may not necessarily reflect the true behavior of the sampled function f_E(z).

	// The integral is calculated as: 
	// < f^2_E(s+\eta_E) > = 1/\sqrt{2\pi}/(\sigma_s^2+\sigma_{in}^2)^{1/2} \times
	// \int_{-\infty}^\infty dx \exp(-x^2/2/(\sigma_s^2+\sigma_{in}^2)) * fe(x)^2

	// Because 0 <= f_E(z) <= 1, we split the integral up into three regions, corresponding
	// to when it is within this range or saturated at the boundaries. The zero boundary does not
	// contribute to the integral, hence we need only compute the other two terms.

	// Get parameters and define convenient names:

	struct interpolantintegrationparams * intparams_fi = (struct interpolantintegrationparams *) params;
	gsl_integration_glfixed_table * t_fi = intparams_fi->GLtable;
	gsl_spline * fi_spline = intparams_fi->splineptr;
	gsl_interp_accel * fi_acc = intparams_fi->accptr;
	int numpts = intparams_fi->numintpts;
	double zi_0 = intparams_fi->zeroroot;
	double zi_1 = intparams_fi->unityroot;
	int ficount = intparams_fi->xcount;
	struct pdf_params modelparams = intparams_fi->modelparams;

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

	// Define some internal variables.

	int fi2_i=0; // loop counter
	double fi2_sum = 0; // Sum to approximate integral.
	double fi2_constpart=0; // Variable to hold constant contribution to integral (from integration region where f_E(z) = 1).
	double weightdummy=0, xdummy=0; // Dummy variables to contain the integration weights and points, respectively.
	
	double zlow = minval(zi_0,zi_1);
	double zhigh = maxval(zi_0,zi_1);

	// Computing the term \int_{z_I1}^\infty dx exp[-x^2/2(\sigma_s^2+\sigma_{in}^2)]/sqrt(2*pi*(sigstim^2+sigin^2)) (closed form expression).
	fi2_constpart = 0.5*(1-gsl_sf_erf(zi_1/sqrt(2*(sigstim*sigstim+sigin*sigin))));

	for(fi2_i=0;fi2_i<numpts*ficount;fi2_i++)
	{
		gsl_integration_glfixed_point (zlow, zhigh, fi2_i, &xdummy, &weightdummy, t_fi);
		// Above command returns integration at (fe2_i)th integration point xdummy, for an integral between
		// the limits of zi_0 and zi_1.

		// Integrate (sum): weight_i*interpolant(x_i)^2*\exp(-0.5*x_i^2/(\sigma_s^2+\sigma_{in}^2))
		double temp_fi = gsl_spline_eval(fi_spline, xdummy, fi_acc);
		fi2_sum += weightdummy*rangethreshold(temp_fi,1.0)*rangethreshold(temp_fi,1.0)*exp(-0.5*xdummy*xdummy/(sigstim*sigstim+sigin*sigin));  
	}
	fi2_sum = fi2_sum/sqrt(2*PI*(sigstim*sigstim+sigin*sigin)); // multiply by constant factors.	
	
	if(zi_0 <= zi_1)
		fi2_sum += fi2_constpart;
	else
		fi2_sum += 1-fi2_constpart;

	return fi2_sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double fegfi_montecarlo_integrand_old(double *k, size_t dim, void *params)
{

	// This is the integral of the two-dimensional integral we need to perform
	// to calculate the average <f_E(s+\eta_E) f_I(-s-\eta_I)>. We have changed variables
	// to x = (\sqrt{2}*u - |b|*w)/2 and y = (\sqrt{2}*u + |b|*w)/2 in order to factorize the
	// function g(x,y) = h(u)h(w). Though the terms f_E(x) and f_I(y) are no longer functions of
	// single variables, this shouldn't present any problems for the integration.

	// Get parameters and define convenient names:
	
	struct fefi_doubleintegralparams * doubleparams = (struct fefi_doubleintegralparams *) params;

	gsl_integration_glfixed_table * t_fi = doubleparams->GLtable_fi;
	gsl_spline * fi_spline = doubleparams->splineptr_fi;
	gsl_interp_accel * fi_acc = doubleparams->accptr_fi;
	double zi_0 = doubleparams->zeroroot_fi;
	double zi_1 = doubleparams->unityroot_fi;
	int ficount = doubleparams->fi_count;

	gsl_integration_glfixed_table * t_fe = doubleparams->GLtable_fe;
	gsl_spline * fe_spline = doubleparams->splineptr_fe;
	gsl_interp_accel * fe_acc = doubleparams->accptr_fe;
	double ze_0 = doubleparams->zeroroot_fe;
	double ze_1 = doubleparams->unityroot_fe;
	int fecount = doubleparams->fe_count;
	
	int numpts = doubleparams->numintpts;
	struct pdf_params modelparams = doubleparams->modelparams;

	//	printf("monte carlo integrand func check: z0 = %g, z1 = %g\n",ze_0,ze_1);

	double stimavg=modelparams.stimavg, sigstim = modelparams.sigstim; 
	double avgin=modelparams.avgin, sigin = modelparams.sigin, vin = modelparams.vin; 
	double avgout=modelparams.avgout, sigout = modelparams.sigout, vout = modelparams.vout;

	double scale = sqrt(2*(sigstim*sigstim+sigin*sigin));
	double nueff = (sigstim*sigstim+sigin*sigin*vin)/(sigstim*sigstim+sigin*sigin);

	double argument_plus = scale*(sqrt((1-nueff)/2)*k[0]+sqrt((1+nueff)/2)*k[1]);
	double argument_minus = scale*(sqrt((1-nueff)/2)*k[0]-sqrt((1+nueff)/2)*k[1]);
	
//	printf("arg+ = %g, arg- = %g\n",argument_plus,argument_minus);
	
	double temp_fe = 0; 
	double temp_fi = 0;
	
	if(ze_0 <= ze_1)
	{
		if(argument_plus < ze_1)
		{
			if(argument_plus > ze_0)
				temp_fe = rangethreshold(gsl_spline_eval(fe_spline, argument_plus, fe_acc),1.0);
			else
				temp_fe = 0;
		}
		else
			temp_fe = 1;
	}
	else
	{
		if(argument_plus < ze_0)
		{
			if(argument_plus > ze_1)
				temp_fe = rangethreshold(gsl_spline_eval(fe_spline, argument_plus, fe_acc),1.0);
			else
				temp_fe = 1;
		}
		else
			temp_fe = 0;
	
	}
		
	//printf("internal checkpoint 1\n");
	
	if(zi_0 <= zi_1)
	{	
		if(argument_minus < zi_1)
		{
			if(argument_minus > zi_0)
				temp_fi = rangethreshold(gsl_spline_eval(fi_spline, argument_minus, fi_acc),1.0);
			else
				temp_fi = 0;
		}
		else
			temp_fi = 1;
	}
	else
	{
		if(argument_minus < zi_0)
		{
			if(argument_minus > zi_1)
				temp_fi = rangethreshold(gsl_spline_eval(fi_spline, argument_minus, fi_acc),1.0);
			else
				temp_fi = 1;
		}
		else
			temp_fi = 0;	
	}
			

	double A = temp_fe*temp_fi*exp(-(k[0]*k[0]+k[1]*k[1]))/PI; 

	return A;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
