function [Hypothesis,Pvalue,KSStat,CriticalValue] = TimeRescalingGoodnessOfFit(EventTimes,Parameters,Delay,FunctionType,Alpha)
% Time rescaling technique to assess goodness of fit of real event times to
% Poisson Process with given parameters
% [Hypothesis,Pvalue,KSStat,CriticalValue] = TimeRescalingGoodnessOfFit(EventTimes,Parameters,Delay,FunctionType,Alpha)
% EventTimes        : real event times (in seconds)
% Parameters        : Poisson Process parameters
% Delay             : Delay parameter of Non Homogeneous Poisson Process
% FunctionType      : intensity function type of Non Homogeneous Poisson Process
% Alpha             : significance level of K-S test for goodness of fit
% Hypothesis        : reject null hypothesis boolean
% Pvalue            : p-value of test
% KSStat            : statistic of KS test, max(abs(F1-F2))
% CriticalValue     : critical value of test (if KSStat>CriticalValue, reject null hypothesis at Alpha significance level)

% compute cumulative intensity function at given points
CumRatePoints = CumulativePoissonRateFunction(Parameters,EventTimes,Delay,FunctionType);

% rescale times (~ exp rvs with mean 1)
RescaledTimes = diff(CumRatePoints);

% transorm times to uniform rvs
TransformedTimes = 1-exp(-RescaledTimes);
% TransformedTimesSorted = sort(TransformedTimes);

% % % RealUniformSamples = ([1:length(TransformedTimesSorted)]-0.5)/length(TransformedTimesSorted);
% % % if Alpha==0.05
% % %     Tolerance=1.36;
% % % elseif Alpha==0.01
% % %     Tolerance=1.63;
% % % end
% % % plot(RealUniformSamples,TransformedTimesSorted);hold on;plot(RealUniformSamples,RealUniformSamples,'r');
% % % plot(RealUniformSamples,RealUniformSamples+Tolerance/sqrt(length(TransformedTimesSorted)),'g');plot(RealUniformSamples,RealUniformSamples-Tolerance/sqrt(length(TransformedTimesSorted)),'g');
% % % hold off;pause;

% plot empirical and known distribution
%figure;cdfplot(TransformedTimes);hold on;plot([0:0.01:1],unifcdf([0:0.01:1],0,1),'r');

% perform K-S test for uniform rv on the interval (0,1)
[Hypothesis Pvalue KSStat CriticalValue]=kstest(TransformedTimes', [TransformedTimes' unifcdf(TransformedTimes',0,1)],Alpha);