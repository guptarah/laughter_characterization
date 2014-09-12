function [LogLikelihood] = ComputeLogLikelihoodPoissonProcess(SamplePoints,RateFunctionPoints)

LogLikelihood = -sum(RateFunctionPoints)+sum(SamplePoints.*log(RateFunctionPoints+eps));%-sum(factorial(SamplePoints));

assert(isreal(LogLikelihood),'Warning: Complex log likelihood result.');