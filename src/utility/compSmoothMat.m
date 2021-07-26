function [smatSum, sMatPred, LapMat] = compSmoothMat(X,Y,isPD)
%% compute the matrix for smoothness penalty
% isPD = 1, output positive define 
% isPD = 0, output semi positive define 
% LapMat : laplacian mat 
%

numSeq = length(X) ; 
sMatPred = cell(1,numSeq) ; 
LapMat = cell(1,numSeq) ; 
smatSum = 0 ; 

PosW = 1 ; 
NegW = 5 ; 

for i =1 : numSeq
    TB = X{i} ; 
    TY = Y(i) ; 
    numFrames = size(TB,2) ; 
    TC = triu(ones(numFrames),-1) - triu(ones(numFrames),2) ; 
    TD = diag(sum(TC,2)) ; 
    TF0 = (TD - TC) ; % for smoothness of beta
    
    if ~isPD
        LapMat{i} = TF0 ; 
        TF = TB * TF0 * TB' ; % for the smoothness of prediction
        sMatPred{i} = TF ;
    else
        % way 1
%         TF1 = Spd_Mat(TF0) ; % Spd_Mat recieves symmetric matrix 
%         TF = TB * TF1 * TB' ; 
%         sMatBeta{i} = TF1 ; 
%         sMatPred{i} = TF ;
        % way 2
        
        TF1 = Spd_Mat(TF0) ;
        TF = TB * TF0 * TB' ; 
        TF = nearestSPD(TF) ; 
        LapMat{i} = TF1 ; 
        sMatPred{i} = TF ;
    end
    
    if TY > 0 
         smatSum = smatSum + TF * PosW  ;
%        smatSum = smatSum + TF * PosW / ( 2 * numFrames)  ;
    else
         smatSum = smatSum + TF * NegW ;
%        smatSum = smatSum + TF * NegW / ( 2 * numFrames) ;
    end
end

eigVal = eig(smatSum) ; 

% if min(eigVal) < 0 || ~isreal(eigVal)
%     error('Hele: Wrong laplacian matrix..., maybe the samples are not enough or the dimension is too large') ; 
% end

