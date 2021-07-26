function [alpha, beta] = monoWeight_down(n,y, isEqual)
%% assign weights to each frame : decrease
% using harmonic function

[alpha, beta] = monoWeight(n,y, isEqual) ; 

beta = beta(end:-1:1) ; 
alpha = alpha(end:-1:1) ; 

