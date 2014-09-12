function [AssignmentIndices] = AssignValuesToHistogramBins(Values,NumOfBins)
% Assigns the values of an input vector to the closest bin of its histogram
% [AssignmentIndices] = AssignValuesToHistogramBins(Values,NumOfBins)
% Values                : vector of input values
% NumOfBins             : number of histogram bins
% AssignmentIndices     : indices of the bin that each value was assigned

if size(Values,1)==1 && size(Values,2)>1
    Values = Values';
end

[n x] = hist(Values,NumOfBins);
[MinVal AssignmentIndices] = min(abs(repmat(Values,[1 length(x)])-repmat(x,[length(Values) 1])),[],2);