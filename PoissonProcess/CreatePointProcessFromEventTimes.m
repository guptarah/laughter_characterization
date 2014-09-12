function [PointProcess EventSequence] = CreatePointProcessFromEventTimes(EventTimes,SamplingRate)
% [PointProcess] = CreatePointProcessFromEventTimes(EventTimes,SamplingRate)
% EventTimes     : times of occuring events (in seconds)
% SamplingRate   : sampling rate
% PointProcess   : time increasing process of events

PointProcess(1:ceil(EventTimes(1)*SamplingRate)-1)=0;
PointProcess(1:ceil(EventTimes(1)*SamplingRate))=1;
for i=2:length(EventTimes);
    StartTime=max(1,ceil(EventTimes(i-1)*SamplingRate));
    EndTime=max(1,ceil(EventTimes(i)*SamplingRate));
    PointProcess(StartTime+1:EndTime)=PointProcess(StartTime)+1;
end

EventSequence=zeros(1,ceil(EventTimes(end)*SamplingRate));
EventSequence(ceil(EventTimes*SamplingRate))=1;