function [alpha, beta] = monoWeight(n,y, isEqual)
%% assign weights to each frame
% using harmonic function

if isEqual
    alpha = ones(n,1) * (1/n) ;
    beta = alpha ;
else
    if y >0 
        lowTri = tril(ones(n)) ; 
        beta0 = 1 / sum(sum(lowTri)) ; 

        beta = beta0 * ones(n,1) ; 
        alpha = lowTri * beta ; 
    else
        alpha = ones(n,1) * (1/n) ;
        beta = alpha ; 
    end
end
