function [GeneratedTimePoints] = GenerateNonHomogeneousPoissonProcessPiecewiseLinear(EventTimes,CumIntensityFunctionEstimator)
% [GeneratedTimePoints] = GenerateNonHomogeneousPoissonProcessPiecewiseLinear(TotalTime,LambdaHomogeneous,LambdaNonHomogeneous,ChangePoints)
% Generates time points of non-homogeneous Poisson Distribution with piecewise linear rate
% EventTimes            : Times (in seconds) of occuring events
% GeneratedTimePoints   : Generated time points of event occurence

TotalRealizations = length(CumIntensityFunctionEstimator);
TotalEvents = length(EventTimes);

i=1;
GeneratedTimePointsHomogeneous = zeros(1,TotalRealizations);
GeneratedTimePoints = zeros(1,TotalRealizations);
GeneratedTimePointsHomogeneous(i) = -log(1-rand);
while GeneratedTimePointsHomogeneous(i)<TotalEvents/TotalRealizations
    m = max(1,floor((TotalEvents+1)*TotalRealizations*GeneratedTimePointsHomogeneous(i)/TotalEvents));
    if m>length(EventTimes)-1
        break
    end
    GeneratedTimePoints(i) = EventTimes(m)+(EventTimes(m+1)-EventTimes(m))*((TotalEvents+1)*TotalRealizations*GeneratedTimePointsHomogeneous(i)/TotalEvents-m);
    i = i+1;
    GeneratedTimePointsHomogeneous(i) = GeneratedTimePointsHomogeneous(i-1)-log(1-rand);
end
m=length(EventTimes)-1;
GeneratedTimePoints(i) = EventTimes(m)+(EventTimes(m+1)-EventTimes(m))*((TotalEvents+1)*TotalRealizations*GeneratedTimePointsHomogeneous(i)/TotalEvents-m);
GeneratedTimePoints(i+1:end)=[];