function [YPoints] = PoissonRateFunction(ParametersToOpt,XPoints,ParametersConstant,FunctionType)
% [YPoints] = PoisonRateFunction(XPoints,Lambda0,Lambda,Alpha,ParametersConstant,FunctionType)
% Computes the rate function of a Poisson process
% YPoints           : Rate function output, 1xN
% XPoints           : points in time axis, 1xN
% ParametersToOpt   : parameter vector
% ParametersConstant: constant parameters of rate function (not for optimization)
% FunctionType      : Type of poisson function
%                     - hom                : homogeneous - constant
%                     - piecewiselin       : piece-wise constant function
%                     - expsum             : sum of time-shifted exponentials
%                     - chisum             : sum of time-shifted chi-square functions
%                     - tukeywinsum        : sum of time-shifted tukey windows
%                     - rectwinsum         : sum of time-shifted rectangle windows
%                     - selfmutual         : self/mutually-exciting Poisson
%                     - selfmutualweighted : self/mutually-exciting Poisson,
%                                            parameters weighted with event magnitude

TukeyWinRatio = 0.5;
MagnHistBins = 10;

if strcmpi(FunctionType,'hom')
    assert(length(ParametersToOpt)==1,'Homogeneous intensity function has only one parameter');
    YPoints = ParametersToOpt*ones(1,length(XPoints));
elseif strcmpi(FunctionType,'piecewiselin')
    % the parameters of the function are the observations of event times
    % and the total number of realizations
    TotalRealizations = ParametersToOpt(end);
    ParametersToOpt(end)=[];
    TotalEvents = length(ParametersToOpt);
    YPoints = zeros(1,length(XPoints));
    for i = 1:length(XPoints)
        PreviousEventIndex = find(ParametersToOpt<=XPoints(i),1,'last');
        if PreviousEventIndex==length(ParametersToOpt)
            YPoints(i) = TotalEvents/((TotalEvents+1)*TotalRealizations*(ParametersToOpt(PreviousEventIndex)-ParametersToOpt(PreviousEventIndex-1)));
        else
            assert(ParametersToOpt(PreviousEventIndex+1)>XPoints(i),'Time after found previous event is still less that time of point to be estimated.');
            YPoints(i) = TotalEvents/((TotalEvents+1)*TotalRealizations*(ParametersToOpt(PreviousEventIndex+1)-ParametersToOpt(PreviousEventIndex)));
        end
    end
elseif strcmpi(FunctionType,'expsum')
    % ParametersToOpt = [Lambda0 Lambda Alpha]
    % Lambda0   : Base Poisson rate parameter
    % Lambda    : Amplitude parameter of exponentials
    % Alpha     : Rate parameter of exponentials
    % YPoints = Lambda0 + Lambda(1)*exp(-Alpha(1)(XPoints-ParametersConstant(1)))*u(Xpoints-ParametersConstant(1)) + ...
    %           Lambda(M)*exp(-Alpha(M)(XPoints-ParametersConstant(M)))*u(Xpoints-ParametersConstant(M))
    assert(mod(length(ParametersToOpt),2)==1,'Number of initial parameters must be odd');
    M = floor(length(ParametersToOpt)/2);
    Lambda0 = ParametersToOpt(1);
    Lambda = ParametersToOpt(2:M+1);
    Alpha = ParametersToOpt(M+2:end);
    YPoints = Lambda0*ones(1,length(XPoints));
    for m=1:M
        IndicesInUse = find(XPoints>ParametersConstant(m));
        YPoints(IndicesInUse) = YPoints(IndicesInUse)+Lambda(m)*exp(-Alpha(m)*(XPoints(IndicesInUse)-ParametersConstant(m)));
    end
elseif strcmpi(FunctionType,'chisum')
    % ParametersToOpt = [Lambda0 Lambda Sigma K]
    % Lambda0   : Base Poisson rate parameter
    % Lambda    : Amplitude parameter of exponentials
    % Sigma     : Width of chi-square
    % K         : Degree of freedom of chi-square
    % YPoints = Lambda0 + Lambda(1)*exp(-Alpha(1)(XPoints-ParametersConstant(1)))*u(Xpoints-ParametersConstant(1)) + ...
    %           Lambda(M)*exp(-Alpha(M)(XPoints-ParametersConstant(M)))*u(Xpoints-ParametersConstant(M))
    assert(mod(length(ParametersToOpt),3)==1,'Number of initial parameters must be 3*M+1');
    M = floor(length(ParametersToOpt)/3);
    Lambda0 = ParametersToOpt(1);
    Lambda = ParametersToOpt(2:M+1);
    Sigma = ParametersToOpt(M+2:2*M+1);
    K = ParametersToOpt(2*M+2:end);
    YPoints = Lambda0*ones(1,length(XPoints));
    for m=1:M
        IndicesInUse = find(XPoints>ParametersConstant(m));
        UnitDrop=chi2pdf(Sigma(m)*(XPoints(IndicesInUse)-ParametersConstant(m)),K(m));
        YPoints(IndicesInUse) = YPoints(IndicesInUse)-Lambda(m)*UnitDrop/max(UnitDrop);
    end
elseif strcmpi(FunctionType,'tukeywinsum')
    % ParametersToOpt = [Lambda0 Lambda Width]
    % Lambda0   : Base Poisson rate parameter
    % Lambda    : Amplitude parameter of tukey windows
    % Width     : Window width (in seconds)
    assert(mod(length(ParametersToOpt),2)==1,'Number of initial parameters must be odd');
    M = floor(length(ParametersToOpt)/2);
    Lambda0 = ParametersToOpt(1);
    Lambda = ParametersToOpt(2:M+1);
    Width = ParametersToOpt(M+2:end);
    YPoints = Lambda0*ones(1,length(XPoints));
    for m=1:M
        IndicesInUse = find(XPoints>ParametersConstant(m));
        if isempty(IndicesInUse)
            continue;
        else
            Window = Lambda(m)*tukeywin(round(Width(m)),TukeyWinRatio)';
            if length(IndicesInUse)>Window
                YPoints(IndicesInUse(1):IndicesInUse(1)+length(Window)-1) = YPoints(IndicesInUse(1):IndicesInUse(1)+length(Window)-1)+Window;
            else
                YPoints(IndicesInUse) = YPoints(IndicesInUse)+Window(1:length(IndicesInUse));
            end
        end
    end
elseif strcmpi(FunctionType,'rectwinsum')
    % ParametersToOpt = [Lambda0 Lambda Width]
    % Lambda0   : Base Poisson rate parameter
    % Lambda    : Amplitude parameter of rectangular windows
    % Width     : Window width (in seconds)
    assert(mod(length(ParametersToOpt),2)==1,'Number of initial parameters must be odd');
    M = floor(length(ParametersToOpt)/2);
    Lambda0 = ParametersToOpt(1);
    Lambda = ParametersToOpt(2:M+1);
    Width = ParametersToOpt(M+2:end);
    YPoints = Lambda0*ones(1,length(XPoints));
    for m=1:M
        IndicesInUse = find(XPoints>ParametersConstant(m));
        if isempty(IndicesInUse)
            continue;
        else
            Window = Lambda(m)*rectwin(round(Width(m)))';
            if length(IndicesInUse)>Window
                YPoints(IndicesInUse(1):IndicesInUse(1)+length(Window)-1) = YPoints(IndicesInUse(1):IndicesInUse(1)+length(Window)-1)+Window;
            else
                YPoints(IndicesInUse) = YPoints(IndicesInUse)+Window(1:length(IndicesInUse));
            end
        end
    end
elseif strcmpi(FunctionType,'selfmutual')
    % ParametersToOpt = [Mu(i) Alpha(1,i) ... Alpha(K,i) Beta(1,i) ... Beta(K,i)], K total processes
    NumOfProcesses = floor(length(ParametersToOpt)/2);
    assert(NumOfProcesses==length(ParametersConstant),'Event times must have as many cells as the number of processes defined by parameters vector');
    Mu = ParametersToOpt(1);
    Alpha = ParametersToOpt(2:NumOfProcesses+1);
    Beta = ParametersToOpt(NumOfProcesses+2:2*NumOfProcesses+1);
    YPoints = Mu*ones(1,length(XPoints));
    for k = 1:NumOfProcesses
        for m=1:length(ParametersConstant{k})
            IndicesInUse = find(XPoints>ParametersConstant{k}(m));
            YPoints(IndicesInUse) = YPoints(IndicesInUse)+Alpha(k)*exp(-Beta(k)*(XPoints(IndicesInUse)-ParametersConstant{k}(m)));
        end
    end
elseif strcmpi(FunctionType,'mutualweighted')
    % ParametersToOpt = [Mu(i) Alpha(1,i) ... Alpha((K-1)*M+K,i) Beta(1,i)... Beta((K-1)*M-K,i)], K total processes, M bins for histogram of events magnitude
    % ParametersConstant = {EventTimes(1) ... EventTimes(K) EventMagnitudes(1) ... EventMagnitudes(K)};
    NumOfProcesses = length(ParametersConstant)/2;
    AssignedIndicesPerProcess = cell(1,NumOfProcesses);
    for k = 1:NumOfProcesses
        assert(length(ParametersConstant{k})==length(ParametersConstant{k+NumOfProcesses}),'Number od SCR event times and magnitudes must be the same')
        AssignedIndicesPerProcess{k} = AssignValuesToHistogramBins(ParametersConstant{NumOfProcesses+k},MagnHistBins);
    end
    Mu = ParametersToOpt(1);
    Alpha = ParametersToOpt(2:(NumOfProcesses-1)*MagnHistBins+NumOfProcesses);
    AlphaSelf = Alpha(1:NumOfProcesses-1); Alpha(1:NumOfProcesses-1)=[];
    Beta = ParametersToOpt((NumOfProcesses-1)*MagnHistBins+NumOfProcesses+1:2*((NumOfProcesses-1)*MagnHistBins+NumOfProcesses-1)+1);
    BetaSelf = Beta(1:NumOfProcesses-1); Beta(1:NumOfProcesses-1)=[];
    YPoints = Mu*ones(1,length(XPoints));
    for k = 1:NumOfProcesses
        for m=1:length(ParametersConstant{k})
            IndicesInUse = find(XPoints>ParametersConstant{k}(m));
            if k==1
                YPoints(IndicesInUse) = YPoints(IndicesInUse)+AlphaSelf(k)*exp(-BetaSelf(k)*(XPoints(IndicesInUse)-ParametersConstant{k}(m)));
            else
                YPoints(IndicesInUse) = YPoints(IndicesInUse)+Alpha(AssignedIndicesPerProcess{k}(m))*exp(-Beta(AssignedIndicesPerProcess{k}(m))*(XPoints(IndicesInUse)-ParametersConstant{k}(m)));
            end
        end
    end
end