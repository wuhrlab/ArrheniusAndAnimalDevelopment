function [numRow,numCol] = calcRectangle(observations)
numRow = floor(sqrt(observations));
numCol = ceil(observations/numRow);
end