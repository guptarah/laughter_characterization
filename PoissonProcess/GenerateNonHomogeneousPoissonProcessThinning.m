function [GeneratedTimePoints] = GenerateNonHomogeneousPoissonProcessThinning(TotalTime,LambdaHomogeneous,LambdaNonHomogeneousParameters,LambdaNonHomogeneousFunction,ChangePoints)
% [GeneratedTimePoints] = GenerateNonHomogeneousPoissonProcessThinning(TotalTime,LambdaHomogeneous,LambdaNonHomogeneous,ChangePoints)
% Generates time points of non-homogeneous Poisson Distribution by thinning
% of a homogeneous Poisson function
% TotalTime                          : Maximum value of time points
% LambdaHomogeneous                  : Rate parameter of homogeneous Poisson process to be thinned
% LambdaNonHomogeneousParameters     : Parameters of nonhomogeneous Poisson process rate
% LambdaNonHomogeneousFunction       : code of nonhomogeneous rate function
% ChangePoints                       : Points where non-homogeneous rate function changes

[GeneratedTimePointsHomogeneous] = GenerateHomogeneousPoissonProcess(LambdaHomogeneous,5*TotalTime);

PoissonRateFunctionValues = PoissonRateFunction(LambdaNonHomogeneousParameters,GeneratedTimePointsHomogeneous,ChangePoints,LambdaNonHomogeneousFunction);

i=1;
k=0;
GeneratedTimePoints = zeros(1,length(GeneratedTimePointsHomogeneous));
while i<=length(GeneratedTimePoints)
    if rand<PoissonRateFunctionValues(i)/LambdaHomogeneous
        k = k+1;
        GeneratedTimePoints(k) = GeneratedTimePointsHomogeneous(i);
    end
    i = i+1;
    if k>0 && k<=length(GeneratedTimePoints) && GeneratedTimePoints(k)>TotalTime
        break
    end
end
if GeneratedTimePoints(k)<TotalTime
    disp('Warning: Time points of non-homogeneous process have not reached desired upper limit.');
end
GeneratedTimePoints(k+1:end) = [];
