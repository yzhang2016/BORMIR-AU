function [predY] = SequenceTest(w,X)
%% 
% augment the feature 

X = [X; ones(1,size(X,2))] ; 
predY = w' * X ; 