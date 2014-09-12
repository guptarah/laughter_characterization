function [Output] = EstimateAndGeneratePoissonProcessPiecewiseLinear(EventsSequence,Config)
% [EstimatedParameters] = EstimateAndGeneratePoissonProcessPiecewiseLinear(EventsSequence,Config)
% EventsSequence    : Sequence of events over time (in samples)
% Config            : Configuration struct containing
%                     - PhysioFS            : sampling frequency
%                     - CellWidth           : duration of each realization (in seconds)
% Output            : Output structure
%                     - GeneratedPointProcess   : simulated point process
%                     - EstimatedParameters     : estimated parameters of Poisson process rate function

assert(mod(length(EventsSequence),Config.PhysioFs)==0,'Sequence of events must have length multiple of the sampling rate');

EventTimes = find(EventsSequence>0)/Config.PhysioFs; % in seconds
PointProcess=cumsum((EventsSequence>0));
EventsPerSecond=sum(reshape(EventsSequence,[Config.PhysioFs length(EventsSequence)/Config.PhysioFs])~=0);
NumberOfEventsPerRealization = sum(buffer(EventsSequence,Config.PhysioFs*Config.CellWidth)~=0);
CumIntensityFunctionEstimator=cumsum(NumberOfEventsPerRealization);

% rate function parameters: [ObservedEventTimes TotalRealizations]
PiecewiseLinFunctionParameters = [EventTimes length(NumberOfEventsPerRealization)];

% generate non-homogeneous Poisson process (event times) with piece-wise linear rate function
[GeneratedTimePoints] = GenerateNonHomogeneousPoissonProcessPiecewiseLinear(EventTimes,CumIntensityFunctionEstimator);
[RateFunctionPoints] = PoissonRateFunction(PiecewiseLinFunctionParameters,[1:length(EventsSequence)/Config.PhysioFs],[],'piecewiselin');

% create point process from time points
GeneratedTimePoints = [0 GeneratedTimePoints];
Slope=Config.CellWidth./[GeneratedTimePoints(1) diff(GeneratedTimePoints)];
Output.GeneratedPointProcess=[];
Level=0;
for i=2:length(GeneratedTimePoints);
    Indices=[max(1,round(GeneratedTimePoints(i-1)*Config.PhysioFs+1)):round(GeneratedTimePoints(i)*Config.PhysioFs)];
    if ~isempty(Indices)
        Output.GeneratedPointProcess(Indices)=Level+Slope(i)*(Indices-Indices(1))/Config.PhysioFs;
        Level=Output.GeneratedPointProcess(Indices(end));
    end
end

if length(Output.GeneratedPointProcess)>length(PointProcess)
    Output.GeneratedPointProcess(length(PointProcess)+1:end)=[];
else
    Output.GeneratedPointProcess(end+1:length(PointProcess))=Output.GeneratedPointProcess(end);
end

% compute RMS error between original and generated Poisson process
Output.RMSError = sqrt(mean((Output.GeneratedPointProcess-PointProcess).^2));
Output.LogLikelihood = ComputeLogLikelihoodPoissonProcess(EventsPerSecond,RateFunctionPoints);
[Output.KSHypothesis Output.KSPvalue Output.KSStat] = TimeRescalingGoodnessOfFit(EventTimes,PiecewiseLinFunctionParameters,[],'piecewiselin',Config.Alpha);