function [Output] = EstimateAndGenerateHomogeneousPoissonProcess(EventsSequence,Config)
% [EstimatedParameters] = EstimateAndGenerateHomogeneousPoissonProcess(EventsSequence,Config)
% EventsSequence    : Sequence of events over time (in samples)
% Config            : Configuration struct containing
%                     - PhysioFS    : sampling frequency
%                     - ZoneLength  : zone length after event occurence in samples
%                     - Alpha       : p-value of hypothesis testing
% Output            : Output structure
%                     - EstimatedParameters     : estimated parameters of Poisson process rate function
%                     - GeneratedPointProcess   : simulated point process
%                     - GeneratedEventsSequence : simulated sequence of events
%                     - RMSError                : root mean square error between original and simulated point process

assert(mod(length(EventsSequence),Config.PhysioFs)==0,'Sequence of events must have length multiple of the sampling rate');

PointProcess=cumsum((EventsSequence>0));
WaitingTimes=diff(find(EventsSequence>0));
EventsPerSecond=sum(reshape(EventsSequence,[Config.PhysioFs length(EventsSequence)/Config.PhysioFs])~=0);
EventTimes = find(EventsSequence>0)/Config.PhysioFs;
%%%%%%%%EventsSequenceZeroPad = EventsSequence(1:floor(length(EventsSequence)/Config.PhysioFs/Config.TimeUnit)*Config.PhysioFs*Config.TimeUnit);
%%%%%%%%EventsPerTimeUnit=sum(reshape(EventsSequenceZeroPad,[Config.PhysioFs*Config.TimeUnit length(EventsSequenceZeroPad)/(Config.PhysioFs*Config.TimeUnit)])~=0);

% empirical CDF of event times
[EmpiricalCdfRealData]=cdfcalc(EventTimes);

% estimate interarrival time of base Poisson Process, that follow
% exponential distributionw with MLE
[Hypothesis Pvalue Output.EstimatedParameters] = TestDataFollowingDistribution(WaitingTimes,'exp',Config.Alpha);
Output.EstimatedParameters = Output.EstimatedParameters/Config.PhysioFs;
RateFunctionPoints = PoissonRateFunction(1/Output.EstimatedParameters,[1:length(EventsPerSecond)],[],'hom');
%%%%%%%%[RateFunctionPoints] = PoissonRateFunction(1/Output.EstimatedParameters,[Config.TimeUnit:Config.TimeUnit:length(EventsPerTimeUnit)*Config.TimeUnit],[],'hom');

% compute likelihood measures
Output.LogLikelihood = ComputeLogLikelihoodPoissonProcess(EventsPerSecond,RateFunctionPoints);
%%%%%%%%Output.LogLikelihood = ComputeLogLikelihoodPoissonProcess(EventsPerTimeUnit,RateFunctionPoints);
Output.AIC = 2*(length(Output.EstimatedParameters)-Output.LogLikelihood);

% compute residual measures
[Output.KSHypothesisME Output.KSPvalueME Output.KSStatME] = KSPlot(EventTimes,1/Output.EstimatedParameters,[],'hom',Config.Alpha);
%[Output.KSHypothesis Output.KSPvalue Output.KSStat] = TimeRescalingGoodnessOfFit(EventTimes,1/Output.EstimatedParameters,[],'hom',Config.Alpha);

% generate random process
Output.RMSError=zeros(1,Config.NumOfIterations);
Output.KSHypothesisSL=zeros(1,Config.NumOfIterations);
Output.KSPvalueSL=zeros(1,Config.NumOfIterations);
Output.KSStatSL=zeros(1,Config.NumOfIterations);
for i=1:Config.NumOfIterations
    % generate homogeneous Poisson process (event times) with estimated parameters
    [GeneratedTimePoints] = GenerateHomogeneousPoissonProcess(1/Output.EstimatedParameters,length(EventsSequence)/Config.PhysioFs);
    % create point process from time points
    [GeneratedPointProcess GeneratedEventsSequence] = CreatePointProcessFromEventTimes(GeneratedTimePoints,Config.PhysioFs);
    if length(GeneratedPointProcess)>length(PointProcess)
        Output.GeneratedPointProcess=GeneratedPointProcess(1:length(PointProcess));
        Output.GeneratedEventsSequence=GeneratedEventsSequence(1:length(PointProcess));
    else
        Output.GeneratedPointProcess=[GeneratedPointProcess GeneratedPointProcess(end)*ones(1,length(PointProcess)-length(GeneratedPointProcess))];
        Output.GeneratedEventsSequence=[GeneratedEventsSequence GeneratedEventsSequence(end)*ones(1,length(PointProcess)-length(GeneratedEventsSequence))];
    end
    % compute RMS error between original and generated Poisson process
    Output.RMSError(i) = sqrt(mean((Output.GeneratedPointProcess-PointProcess).^2));
    % compute empirical CDF from simulated data
    [EmpiricalCdfSimulatedData]=cdfcalc(GeneratedTimePoints);
    % perform KS-test between CDFs
    [Output.KSHypothesisSL(i) Output.KSPvalueSL(i) Output.KSStatSL(i)] = kstest2(EmpiricalCdfRealData,EmpiricalCdfSimulatedData,Config.Alpha);
end
