#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

struct paramsets
{
	double sigstim;
	double sigin;
	double sigout;
	double vin;
	double vout;
};

int main(void)
{

	// Define parameter sets
	
	double sig_s_params[] = {1.0};
	double sig_in_params[] = {2.0}; //{0.1,1.0};
	double sig_out_params[] = {0.1,0.25,0.5,0.75,1.0};
	double vin_params[] = {-0.25, 0.0625, 0.3750, 0.6875, 1.0}; //gives v_eff = {0,0.25,0.5,0.75,1.0};
	double vout_params[] = {-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0}; //{-0.95,-0.5,-0.1,0.0,0.1,0.5,0.95};
	double kappa_params[] = {0.0,0.25,0.5,0.75,1.0}; 
	
	int seedvalues[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946}; //, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269};
	
	int OOvOF[] = {0,1};
	// Get number of elements in each parameter set.
	
	int num_sig_s = (int)(sizeof(sig_s_params)/sizeof(double));
	int num_sig_in = (int)(sizeof(sig_in_params)/sizeof(double));
	int num_sig_out = (int)(sizeof(sig_out_params)/sizeof(double));
	int num_vin = (int)(sizeof(vin_params)/sizeof(double));
	int num_vout = (int)(sizeof(vout_params)/sizeof(double));
	int num_kappa = (int)(sizeof(kappa_params)/sizeof(double));
	int num_seed = (int)(sizeof(seedvalues)/sizeof(int));
	int num_OOvOF = (int)(sizeof(OOvOF)/sizeof(int));
	
	// Set file pointers and names, then open output file
	
	FILE *output;
	char outputname[100];
	sprintf(outputname,"%s","sweepfile");

	output = fopen(outputname,"w");
	
	// loop indicies
	int ssi=0, sii = 0, soi = 0, vii = 0, voi = 0, ki = 0, sdi = 0, oovofi = 0; 

	for(oovofi=0;oovofi < num_OOvOF; oovofi++)
	{
		for(ssi=0;ssi<num_sig_s;ssi++)
		{
			for(sii=0;sii<num_sig_in;sii++)
			{
				for(soi=0;soi<num_sig_out;soi++)
				{
					for(vii=0; vii < num_vin; vii++)
					{
						for(voi=0; voi < num_vout; voi++)
						{
							for(ki=0; ki < num_kappa; ki++)
							{
								for(sdi=0;sdi<num_seed;sdi++)
								{
									//fprintf(output,"%s %g %g %g %g %g %g %i\n","LD_LIBRARY_PATH=/gscratch/riekesheabrown/local/lib/ /gscratch/riekesheabrown/bradenb/ppie_solver3/executables/ppie_solver3",sig_s_params[ssi],sig_in_params[sii],sig_out_params[soi],vin_params[vii],vout_params[voi],kappa_params[ki],seedvalues[sdi]);
									fprintf(output,"%s %g %g %g %g %g %g %i %i\n","/gscratch/riekesheabrown/bradenb/ppie_solver4/executables/ppie_solver4",sig_s_params[ssi],sig_in_params[sii],sig_out_params[soi],vin_params[vii],vout_params[voi],kappa_params[ki],seedvalues[sdi],OOvOF[oovofi]);

								}
							}
						}
					}
				}
			}
		}
	}
	
	// End of program clean up
	fclose(output);
	
	return 0;
}
