function [Hypothesis Pvalue DistributionParameters] = TestDataFollowingDistribution(Data,Distribution,Alpha)

if size(Data,1)>1 && size(Data,2)>1
    error('Data must be a vector');
end

if size(Data,1)==1 && size(Data,2)>1
    Data = Data';
end

if nargin==2 || isempty(Alpha)
    Alpha = 0.1;
end

if strcmp(Distribution,'gamma')
    [DistributionParameters] = gamfit(Data);
    [Hypothesis,Pvalue]=kstest(Data, [Data gamcdf(Data,DistributionParameters(1),DistributionParameters(2))],Alpha);
elseif strcmp(Distribution,'exp')
    [DistributionParameters] = expfit(Data);
    [Hypothesis,Pvalue]=kstest(Data, [Data expcdf(Data,DistributionParameters)],Alpha);
elseif strcmp(Distribution,'gp')
    [DistributionParameters] = gpfit(Data);
    if sum(isnan(DistributionParameters))>0
        Hypothesis=0;
        Pvalue=1;
    else
        [Hypothesis,Pvalue]=kstest(Data, [Data gpcdf(Data,DistributionParameters(1),DistributionParameters(2),0)],Alpha);
    end
elseif strcmp(Distribution,'normal')
    [DistributionParameters(1) DistributionParameters(2)] = normfit(Data);
    [Hypothesis,Pvalue]=kstest(Data, [Data normcdf(Data,DistributionParameters(1),DistributionParameters(2))],Alpha);
% elseif strcmp(Distribution,'binomial')
%     [DistributionParameters] = binofit(sum(Data==1),length(Data));
%     [Hypothesis,Pvalue]=kstest(Data, [Data binocdf([1:length(Data)]',length(Data),DistributionParameters)],Alpha);
elseif strcmp(Distribution,'poisson')
    [DistributionParameters] = poissfit(Data,Alpha);
    [Hypothesis,Pvalue]=kstest(Data, [Data poisscdf(Data,DistributionParameters)],Alpha);
else
    error('Given distribution not supported yet.');
end