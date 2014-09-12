The main functions of the code are the following. They all contains comments about their functionality, input and output.

- EstimateAndGeneratePoissonProcess.m		: Estimates the parameters of a non-homogeneous Poisson Process (PP) and generates samples following this process.
- EstimateAndGeneratePoissonProcess2.m		: Estimates the parameters of a non-homogeneous Poisson Process (PP) and generates samples following this process. If the interval length of parameter estimation is not the unit time (e.g. 1sec, like it was in my paper), use this function to estimate the PP parameters.
- EstimateAndGenerateHomogeneousPoissonProcess.m	: Estimates the parameter of a homogeneous PP and generates samples following this process.
- EstimateAndGeneratePoissonProcess2d.m		: Estimates the parameters of a 2-D PP. I have started writing it, but because I didn't need it at the end, it is not complete. Tell me if you need it, so I can finish it.


The following functions are called from the main functions that were described above. Here is a brief explanation for the functionality of each of them.
- GenerateHomogeneousPoissonProcess.m		: Generates a homogeneous PP.
- GenerateNonHomogeneousPoissonProcessThinning.m: Generates a non-homogeneous PP with the method of thinning (Lewis and Schedler 1979).
- GenerateNonHomogeneousPoissonProcessPiecewiseLinear.m : Generates a piece-wise linear PP.
- InitializePoissonRateFunctionParameters.m	: Initializes the PP rate function parameters with appropriate values for the JASPER data. You might need to change these values.
- PoissonRateFunction.m				: Computes different kinds of PP rate functions. Not sure if every kind of function works, because I updated the code many times. The 'expsum' preference (sum of exponentials) works for sure.
- CumulativePoissonRateFunction.m		: Computes the cumulative rate function of the PP (as above, it works for the expsum choice for sure).
- PoissonRateFunction2d.m			: 2D PP rate function. Not complete yet.
- ComputeLogLikelihoodPoissonProcess.m		: Computes the log-likelihood of the PP.
- KSPlot.m					: KS plot statistics to assess goodness of fit.
- TimeRescalingGoodnessOfFit.m			: Goodness of fit measures with Time Rescaling (Brown et. al 2001, not the exact reference, but will give you an idea)

Run ToyExample.m to as a first example. This uses SCROccurencesExample.mat and ChangePoints.mat files.
