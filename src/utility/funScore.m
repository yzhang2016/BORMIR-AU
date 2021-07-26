function [score] = funScore(w,alpha,beta,LapMat,Y,X,FMat,lambda,gamma,rho,isWeighted,normP)
%% compute the objective function score 
% LapMat : laplacian mat
%
OX = cell(size(alpha)) ; % average instance
for i = 1 : length(OX)
    OX{i} = X{i}*alpha{i} ; 
end
[OX,numVec] = cellToMat(OX) ; 
OY = Y(:) ;

sampW = ones(length(OY),1) ;
if isWeighted
    numNeu = sum(OY == 0); 
    numPos = length(OY) - numNeu ; 
%     sampW(OY == 0) = sqrt(1 / numNeu) ; 
%     sampW(OY ~= 0) = sqrt(1 / numPos) ; 
    sampW(OY == 0) = 1 ; 
    sampW(OY ~= 0) = sqrt(numNeu / numPos) ; 
end

[OBeta,~] = cellToMat(beta') ; 
OD = diag(sampW.*sampW) ;

A = (OY - OX' * w)' * OD * (OY - OX' * w) ; 
B =  w' * w ; 
C = w' * FMat * w ; 

DA = 0 ; 
DB = 0 ; 
for i =1 : length(LapMat)
    TF = LapMat{i} ;
    TA = alpha{i} ; 
    TB = beta{i} ; 
    DA = DA + TA' * TF * TA ; 
    DB = DB + TB' * TF * TB ;
end

% D =  sum(abs(OBeta).^normP) ;

score = A + lambda * B  + rho * C + gamma * DA; 

