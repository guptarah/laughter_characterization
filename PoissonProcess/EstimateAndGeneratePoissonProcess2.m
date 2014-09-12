function [Output] = EstimateAndGeneratePoissonProcess2(EventsSequence,Config)
% [Output] = EstimateAndGeneratePoissonProcess(EventsSequence,Config)
% EventsSequence    : Sequence of events over time (in samples)
% Config            : Configuration struct containing
%                     - PhysioFS           : sampling frequency
%                     - RateFunctionName   : code name of rate function (e.g. expsum)
%                     - InitialParameters  : initial rate function parameters and low/high boundaries (3xNumOfParameters)
%                     - ChangePoints       : points (in samples) where rate function changes
%                     - NumOfIterations    : number of iterations for generating random process
%                     - TimeUnit           : interval length of parameter estimation (in seconds)
% Output            : Output structure
%                     - GeneratedPointProcess   : simulated point process (last iteration)
%                     - GeneratedEventsSequence : simulated sequence of events (last iteration)
%                     - EstimatedParameters     : estimated parameters of Poisson process rate function
%                     - RMSError                : root mean square error between original and simulated point process (vector with all iterations)
%                     - KSStat                  : statistic of KS goodness of fit test
%                     - KSPvalue                : p-value of KS goodness of fit test
%                     - LogLikelihood           : log likelihood of poisson process with estimated parameters
%                     - AIC                     : AIC of poisson process with estimated parameters

assert(mod(length(EventsSequence),Config.PhysioFs)==0,'Sequence of events must have length multiple of the sampling rate');

% Config.ChangePoints=Config.ChangePoints/TimeUnit;
PointProcess=cumsum((EventsSequence>0));
EventsSequenceZeroPad = EventsSequence(1:floor(length(EventsSequence)/Config.PhysioFs/Config.TimeUnit)*Config.PhysioFs*Config.TimeUnit);
EventsPerTimeUnit=sum(reshape(EventsSequenceZeroPad,[Config.PhysioFs*Config.TimeUnit length(EventsSequenceZeroPad)/(Config.PhysioFs*Config.TimeUnit)])~=0);
EventTimes = find(EventsSequence>0)/Config.PhysioFs;

% empirical CDF of event times
[EmpiricalCdfRealData]=cdfcalc(EventTimes);

% estimate rate parameters with least squares fit
f = @(x,y)PoissonRateFunctionNew(x,y,Config.ChangePoints,Config.RateFunctionName);
[Output.EstimatedParameters]=lsqcurvefit(f,Config.InitialParameters(1,:),[1:length(EventsPerTimeUnit)],EventsPerTimeUnit);
%[Output.EstimatedParameters]=lsqcurvefit(f,Config.InitialParameters(1,:),[1:length(EventsPerTimeUnit)]*Config.TimeUnit,EventsPerTimeUnit/((length(EventsSequence)/Config.PhysioFs)/length(EventsPerTimeUnit)),Config.InitialParameters(2,:),Config.InitialParameters(3,:));
[RateFunctionPoints] = PoissonRateFunctionNew(Output.EstimatedParameters,[Config.TimeUnit:Config.TimeUnit:length(EventsPerTimeUnit)*Config.TimeUnit],Config.ChangePoints,Config.RateFunctionName);
RateFunctionPoints(RateFunctionPoints<0) = eps;
Output.RateFunctionPoints = RateFunctionPoints;
% compute likelihood measures
Output.LogLikelihood = ComputeLogLikelihoodPoissonProcess(EventsPerTimeUnit,RateFunctionPoints);
Output.AIC = 2*(length(Output.EstimatedParameters)-Output.LogLikelihood);

% compute residual measures
[Output.KSHypothesisME Output.KSPvalueME Output.KSStatME] = KSPlot(EventTimes,Output.EstimatedParameters,Config.ChangePoints,Config.RateFunctionName,Config.Alpha);
%[Output.KSHypothesis Output.KSPvalue Output.KSStat] = TimeRescalingGoodnessOfFit(EventTimes,Output.EstimatedParameters,Config.ChangePoints,Config.RateFunctionName,Config.Alpha);

%% generate random process
%Output.RMSError=zeros(1,Config.NumOfIterations);
%Output.KSHypothesisSL=zeros(1,Config.NumOfIterations);
%Output.KSPvalueSL=zeros(1,Config.NumOfIterations);
%Output.KSStatSL=zeros(1,Config.NumOfIterations);
%for i=1:Config.NumOfIterations
%    % generate non-homogeneous Poisson process (event times) with estimated parameters
%    [GeneratedTimePoints] = GenerateNonHomogeneousPoissonProcessThinning(length(EventTimes),max(RateFunctionPoints),Output.EstimatedParameters,Config.RateFunctionName,Config.ChangePoints);
%    % create point process from time points
%    [GeneratedPointProcess GeneratedEventsSequence] = CreatePointProcessFromEventTimes(GeneratedTimePoints,Config.PhysioFs);
%    if length(GeneratedPointProcess)>length(PointProcess)
%        Output.GeneratedPointProcess=GeneratedPointProcess(1:length(PointProcess));
%        Output.GeneratedEventsSequence=GeneratedEventsSequence(1:length(PointProcess));
%    else
%        Output.GeneratedPointProcess=[GeneratedPointProcess GeneratedPointProcess(end)*ones(1,length(PointProcess)-length(GeneratedPointProcess))];
%        Output.GeneratedEventsSequence=[GeneratedEventsSequence GeneratedEventsSequence(end)*ones(1,length(PointProcess)-length(GeneratedEventsSequence))];
%    end
%    % compute RMS error between original and generated Poisson process
%    Output.RMSError(i) = sqrt(mean((Output.GeneratedPointProcess-PointProcess).^2));
%    % compute empirical CDF from simulated data
%    [EmpiricalCdfSimulatedData]=cdfcalc(GeneratedTimePoints);
%    % perform KS-test between CDFs
%    [Output.KSHypothesisSL(i) Output.KSPvalueSL(i) Output.KSStatSL(i)] = kstest2(EmpiricalCdfRealData,EmpiricalCdfSimulatedData,Config.Alpha);
%end
