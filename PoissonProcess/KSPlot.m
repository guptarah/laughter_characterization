function [Hypothesis,Pvalue,KSStat] = KSPlot(EventTimes,Parameters,Delay,FunctionType,Alpha)
% KS plot to assess goodness of fit of real event times to
% Poisson Process with given parameters
% [Hypothesis,Pvalue,KSStat,CriticalValue] = KSPlot(EventTimes,Parameters,Delay,FunctionType,Alpha)
% EventTimes        : real event times (in seconds)
% Parameters        : Poisson Process parameters
% Delay             : Delay parameter of Non Homogeneous Poisson Process
% FunctionType      : intensity function type of Non Homogeneous Poisson Process
% Alpha             : significance level of K-S test for goodness of fit
% Hypothesis        : reject null hypothesis boolean
% Pvalue            : p-value of test
% KSStat            : statistic of KS test, max(abs(F1-F2))
% CriticalValue     : critical value of test (if KSStat>CriticalValue, reject null hypothesis at Alpha significance level)

% calculate empirical CDF
[ycdf xcdf]=cdfcalc(EventTimes);
EmpiricalCdf=ycdf(1:end-1);

% compute model CDF
CumRatePoints = CumulativePoissonRateFunction(Parameters,xcdf,Delay,FunctionType);
ModelCdf=CumRatePoints/max(CumRatePoints);

%figure;plot(EmpiricalCdf,ModelCdf); % must be 45 degree line

% perform K-S test between empirical and model CDF
% Hypothesis=0 : the two sample vectors come from the same distribution
[Hypothesis Pvalue KSStat]=kstest2(EmpiricalCdf,ModelCdf,Alpha);