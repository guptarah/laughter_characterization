function [ZPoints] = PoissonRateFunction2d(Parameters,XYPoints,Delay,FunctionType)
% [ZPoints] = PoisonRateFunction(XYPoints,Lambda0,Lambda,Alpha,Delay,FunctionType)
% Computes the rate function of a Poisson process
% ZPoints               : Rate function output, 1xN
% XYPoints              : points in time-space axis, 1xN
% Parameters            : parameter vector 
% Delay                 : Delay parameter of exponentials in space and time, 1xM
% FunctionType          : Type of poisson function
%                         - expsum : sum of time-shifted exponentials

assert(sum(NumOfPointsInXYAxis)==length(XYPoints),'Defined and actual number of points in x and y axis are not equal');

if strcmpi(FunctionType,'expsum')
    % Parameters = [Lambda0 Lambda Alpha Beta]
    % Lambda0   : Base Poisson rate parameter
    % Lambda    : Amplitude parameter of exponentials
    % Alpha     : Rate parameter of exponentials in time
    % Beta      : Rate parameter of exponentials in space
    assert(mod(length(Parameters),3)==1,'Number of initial parameters must 3*M+1');
    M = floor(length(Parameters)/3);
    assert(size(Delay,1)==2 && size(Delay,2)==M,'Number of parameters and offset times must be equal');
    Lambda0 = Parameters(1);
    Lambda = Parameters(2:M+1);
    Alpha = Parameters(M+2:2*M+1);
    Beta = Parameters(2*M+2:end);
    ZPoints = Lambda0*ones(1,size(XYPoints,2));
    for m=1:M
        IndicesInUse = intersect(find(XYPoints(1,:)>Delay(1,m)),find(XYPoints(2,:)>Delay(2,m)));
        ZPoints(IndicesInUse) = ZPoints(IndicesInUse)+Lambda(m)*exp(-Alpha(m)*(XYPoints(1,IndicesInUse)-Delay(1,m))-Beta(m)*(XYPoints(2,IndicesInUse)-Delay(2,m)));
    end
end