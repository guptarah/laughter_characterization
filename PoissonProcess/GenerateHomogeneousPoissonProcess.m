function [GeneratedTimePoints] = GenerateHomogeneousPoissonProcess(Lambda,TotalTime)
% [GeneratedTimePoints] = GenerateHomogeneousPoissonProcess(Lambda,TotalTime)
% generate points of homegeneous Poisson process with parameter Lambda in
% the interval [0,TotalTime]
% Lambda                : Rate parameter of homogeneous Poisson process
% TotalTime             : Maximum time of generated points
% GeneratedTimePoints   : Generated points

TimeOfNewEvent = -log(rand)/Lambda;
GeneratedTimePoints = [];
while TimeOfNewEvent<TotalTime
    GeneratedTimePoints(end+1) = TimeOfNewEvent;
    TimeOfNewEvent = TimeOfNewEvent-log(rand)/Lambda;
end