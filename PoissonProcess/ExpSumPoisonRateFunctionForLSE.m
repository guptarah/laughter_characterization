function [YPoints] = ExpSumPoisonRateFunctionForLSE(Parameters,XPoints)
% YPoints           : Rate function output, 1xN
% Parameters        : [l(0) l(1) .... l(M) a(1) a(2) ... a(M)]
% XPoints           : points in time axis, 1xN

global ChangePoints;

assert(mod(length(Parameters),2)==1,'Number of initial parameters must be odd');

[YPoints] = ExpSumPoisonRateFunction(Parameters,XPoints,ChangePoints);