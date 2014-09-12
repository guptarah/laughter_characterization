function [Output] = EstimateAndGeneratePoissonProcess2d(EventsSequence,XYPoints,Config)
% [EstimatedParameters] = EstimateAndGeneratePoissonProcess2d(EventsSequence,Config)
% EventsSequence    : Sequence of events over time (in samples)
% XYPoints          : 
% Config            : Configuration struct containing
%                     - PhysioFS           : sampling frequency
%                     - RateFunctionName   : code name of rate function (e.g. expsum)
%                     - InitialParameters  : initial rate function parameters and low/high boundaries (3xNumOfParameters)
%                     - ChangePoints       : points (in samples) where rate function changes
% Output            : Output structure
%                     - GeneratedPointProcess   : simulated point process
%                     - GeneratedEventsSequence : simulated sequence of events
%                     - RMSError                : root mean square error between original and simulated point process
%                     - EstimatedParameters     : estimated parameters of Poisson process rate function

assert(mod(length(EventsSequence),Config.PhysioFs)==0,'Sequence of events must have length multiple of the sampling rate');

PointProcess=cumsum((EventsSequence>0));
EventsPerSecond=sum(reshape(EventsSequence,[Config.PhysioFs length(EventsSequence)/Config.PhysioFs])~=0);
EventTimes = find(EventsSequence>0)/Config.PhysioFs;

% estimate rate parameters with least squares fit
f = @(x,y)PoissonRateFunction(x,y,Config.ChangePoints,Config.RateFunctionName);
%[Output.EstimatedParameters]=lsqcurvefit(f,Config.InitialParameters(1,:),[1:length(EventsPerSecond)],EventsPerSecond);
[Output.EstimatedParameters]=lsqcurvefit(f,Config.InitialParameters(1,:),[1:length(EventsPerSecond)],EventsPerSecond,Config.InitialParameters(2,:),Config.InitialParameters(3,:));
[RateFunctionPoints] = PoissonRateFunction(Output.EstimatedParameters,[1:length(EventsPerSecond)],Config.ChangePoints,Config.RateFunctionName);

% generate non-homogeneous Poisson process (event times) with estimated parameters
[GeneratedTimePoints] = GenerateNonHomogeneousPoissonProcessThinning(length(EventsPerSecond),max(RateFunctionPoints),Output.EstimatedParameters,Config.RateFunctionName,Config.ChangePoints);

% create point process from time points
[Output.GeneratedPointProcess Output.GeneratedEventsSequence] = CreatePointProcessFromEventTimes(GeneratedTimePoints,Config.PhysioFs);

if length(Output.GeneratedPointProcess)>length(PointProcess)
    Output.GeneratedPointProcess(length(PointProcess)+1:end)=[];
    Output.GeneratedEventsSequence(length(PointProcess)+1:end)=[];
else
    Output.GeneratedPointProcess(end+1:length(PointProcess))=Output.GeneratedPointProcess(end);
    Output.GeneratedEventsSequence(end+1:length(PointProcess))=Output.GeneratedEventsSequence(end);
end


% compute RMS error between original and generated Poisson process
Output.RMSError = sqrt(mean((Output.GeneratedPointProcess-PointProcess).^2));
Output.LogLikelihood = ComputeLogLikelihoodPoissonProcess(EventsPerSecond,RateFunctionPoints);
[Output.KSHypothesis Output.KSPvalue Output.KSStat] = TimeRescalingGoodnessOfFit(EventTimes,Output.EstimatedParameters,Config.ChangePoints,Config.RateFunctionName,Config.Alpha);