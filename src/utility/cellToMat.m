function [outMat,numVec] = cellToMat(inCell)
%% covert mat to cell 
outMat = cell2mat(inCell) ; 
for i =1 : length(inCell)
    numVec(i) = size(inCell{i},2) ; 
end