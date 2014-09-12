function [YPoints] = CumulativePoissonRateFunction(Parameters,XPoints,Delay,FunctionType)
% [YPoints] = CumulativePoissonRateFunction(XPoints,Lambda0,Lambda,Alpha,Delay,FunctionType)
% Computes the cumulative rate function of a Poisson process
% YPoints           : Rate function output, 1xN
% XPoints           : points in time axis, 1xN
% Parameters        : parameter vector
% Delay             : Delay parameter of exponentials
% FunctionType      : Type of poisson function
%                     - hom         : homogeneous - constant
%                     - piecewiselin: pieciewise constant function
%                     - expsum      : sum of time-shifted exponentials
%                     - rectwinsum  : sum of time-shifted rectangle windows

if strcmpi(FunctionType,'hom')
    assert(length(Parameters)==1,'Homogeneous intensity function has only one parameter');
    YPoints = Parameters*XPoints;
elseif strcmpi(FunctionType,'piecewiselin')
    % the parameters of the function are the observations of event times
    % and the total number of realizations
    TotalRealizations = Parameters(end);
    Parameters(end)=[];
    TotalEvents = length(Parameters);
    YPoints = zeros(1,length(XPoints));
    for i = 1:length(XPoints)
        NextEventIndex = find(Parameters>=XPoints(i),1,'first');
        if NextEventIndex==1
            YPoints(i) = NextEventIndex*TotalEvents/((TotalEvents+1)*TotalRealizations)+TotalEvents*(XPoints(i)-Parameters(NextEventIndex))/((TotalEvents+1)*TotalRealizations*(Parameters(NextEventIndex+1)-Parameters(NextEventIndex)));
        else
            assert(Parameters(NextEventIndex-1)<XPoints(i),'Time after found previous event is still less that time of point to be estimated.');
            YPoints(i) = (NextEventIndex-1)*TotalEvents/((TotalEvents+1)*TotalRealizations)+TotalEvents*(XPoints(i)-Parameters(NextEventIndex-1))/((TotalEvents+1)*TotalRealizations*(Parameters(NextEventIndex)-Parameters(NextEventIndex-1)));
        end
    end
elseif strcmpi(FunctionType,'expsum')
    % Parameters = [Lambda0 Lambda Alpha]
    % Lambda0   : Base Poisson rate parameter
    % Lambda    : Amplitude parameter of exponentials
    % Alpha     : Rate parameter of exponentials
    assert(mod(length(Parameters),2)==1,'Number of initial parameters must be odd');
    M = floor(length(Parameters)/2);
    Lambda0 = Parameters(1);
    Lambda = Parameters(2:M+1);
    Alpha = Parameters(M+2:end);
    YPoints = Lambda0*XPoints;
    for m=1:M
        IndicesInUse = find(XPoints>Delay(m));
        YPoints(IndicesInUse) = YPoints(IndicesInUse)+Lambda(m)/Alpha(m)*(1-exp(-Alpha(m)*(XPoints(IndicesInUse)-Delay(m))));
    end
elseif strcmpi(FunctionType,'rectwinsum')
    % Parameters = [Lambda0 Lambda Width]
    % Lambda0   : Base Poisson rate parameter
    % Lambda    : Amplitude parameter of rectangular windows
    % Width     : Window width (in seconds)
    assert(mod(length(Parameters),2)==1,'Number of initial parameters must be odd');
    M = floor(length(Parameters)/2);
    Lambda0 = Parameters(1);
    Lambda = Parameters(2:M+1);
    Width = Parameters(M+2:end);
    YPoints = Lambda0*XPoints;
    for m=1:M
        IndicesInUse1 = find(XPoints>Delay(m));
        IndicesInUse2 = find(XPoints>Delay(m)+Width(m));
        YPoints(IndicesInUse1) = YPoints(IndicesInUse1)+Lambda(m)*(XPoints(IndicesInUse1)-Delay(m));
        YPoints(IndicesInUse2) = YPoints(IndicesInUse2)-Lambda(m)*(XPoints(IndicesInUse2)-Delay(m)-Width(m));
    end
elseif strcmpi(FunctionType,'selfmutual')
    % [Mu(i) Alpha(1,i) ... Alpha(K,i) Beta(1,i) ... Beta(K,i)], K total processes
    NumOfProcesses = floor(length(Parameters)/2);
    assert(NumOfProcesses==length(Delay),'Event times must have as many cells as the number of processes defined by parameters vector');
    Mu = Parameters(1);
    Alpha = Parameters(2:NumOfProcesses+1);
    Beta = Parameters(NumOfProcesses+2:2*NumOfProcesses+1);
    YPoints = Mu*XPoints;
    for k = 1:NumOfProcesses
        for m=1:length(Delay{k})
            IndicesInUse = find(XPoints>Delay{k}(m));
            YPoints(IndicesInUse) = YPoints(IndicesInUse)+(Alpha(k)/Beta(k))*(1-exp(-Beta(k)*(XPoints(IndicesInUse)-Delay{k}(m))));
        end
    end
end