function [outCell] = matToCell(inmat,numVec)
%% covert mat to cell 
% each column is an observation
outCell = mat2cell(inmat,1,numVec) ; 