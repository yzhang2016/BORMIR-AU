function [w,alpha,beta,score,WRes,alphaRes,alphaResS] = BOMIR_IV_B_ext(Y,Y0,X,lambda,gamma,rho,yita,augRate,option)
%% constrained multi-instance regression : two beta
%
% two beta, one for positive order and the other for negative order
% solve the two beta toghther  
%  
%

if isempty(option)
    maxIter = 30 ; 
    thresh = 1e-4 ;
    isWeighted = 0 ; 
    isEqualFrameW = 0 ; 
else
    maxIter = option.maxIter ; 
    thresh = option.thresh ; 
    isWeighted = option.isWeighted ;
    isEqualFrameW = option.isEqualFrameW ;
    isRmvOrdinal = option.isRmvOrdinal ;
end

numSeq = length(Y) ; 

% augment features 
alpha = cell(1,numSeq) ; 
beta = cell(1,numSeq) ; 

for i = 1 : numSeq
    X{i} = [X{i}; ones(1,size(X{i},2))] ; 
    numFrame = size(X{i},2) ; 
    [temA, temB] = monoWeight(numFrame,Y(i),isEqualFrameW) ; 
    [temAS,temBS] = monoWeight_down(numFrame,Y(i),isEqualFrameW) ; 
    alpha{i} = temA ; 
    beta{i} = temB ; 
    alphaS{i} = temAS ; 
    betaS{i} = temBS ;
end

isPD = 0 ; 
[FMat, ~ ,LapMat] = compSmoothMat(X,Y,isPD); 

n = size(X{i},1) ; % augmented dim

WRes = [] ; 
alphaRes = {alpha} ; 
alphaResS = {alphaS} ; 

w = zeros(n,1) ; 
score = [] ;  

% matlab quadratic programming
options = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','off');
    
temScore = funScore(w,alpha,alphaS,beta,betaS,LapMat,Y,Y0,X,FMat,...
                        lambda,gamma,rho,isWeighted,augRate) ; 
score = [score;temScore] ; 

for TT = 1 : maxIter
    
    % Subproblem: (1) Fix b_i, optimize w
    OX = cell(size(alpha)) ; % average instance
    OXS = cell(size(alphaS)) ; % average instance
    for i = 1 : length(OX)
        OX{i} = X{i}*alpha{i} ; 
        OXS{i} = X{i}*alphaS{i} ; 
    end
    [OX,numVec] = cellToMat(OX) ; 
    [OXS,numVec] = cellToMat(OXS) ;
    OY = Y(:) ;
    OYS = Y0(:) ; 

    sampW = ones(length(OY),1) ;
    if isWeighted
        numNeu = sum(OY == 0); 
        numPos = length(OY) - numNeu ; 
    %     sampW(OY == 0) = sqrt(1 / numNeu) ; 
    %     sampW(OY ~= 0) = sqrt(1 / numPos) ; 
        sampW(OY == 0) = 1 ; 
        sampW(OY ~= 0) = sqrt(numNeu / numPos) ; 
    end

    
    % ===============================================
    % Matlab quadprog 
%     augRate = 1; 
    
    epsMat = eye(size(OX,1))* 1e-12 ; 
    OD = diag(sampW.*sampW) ;
    OH = 2 * (OX*OD*OX' + augRate* OXS*OD*OXS' + lambda * eye(n) + rho * FMat+ epsMat) ; 
    Of = -2 * (OX*OD*OY + augRate* OXS*OYS) ; 
    OW0 = w ; 
    OW = quadprog(OH,Of,[],[],[],[],[],[],OW0,options) ; 
    % ===============================================
    
    w = OW ; 

    temScore = funScore(w,alpha,alphaS,beta,betaS,LapMat,Y,Y0,X,FMat,...
                            lambda,gamma,rho,isWeighted,augRate) ; 
    score = [score;temScore] ; 
    
    % Subproblem: (2)  Fix w, optimize b_i 
    parfor i = 1 : numSeq
        TY = Y(i) ; 
        TY0 = Y0(i); 
        TX = X{i} ;
        numFrame = size(TX,2) ;    

        if isRmvOrdinal
            TL = eye(numFrame) ;
        else
            if TY > 0
                TL = tril(ones(numFrame)) ;
            else
                TL = eye(numFrame) ;
            end 
        end
        % ===============================================
        TF = LapMat{i};  
        
         TF = TL' * TF * TL ;  % add smoothness on alpha (the salience)
         TF2 = TL * TF * TL' ;  % add smoothness on alpha (the salience)
         
%        TF = TF ; % add smoothness on beta (the change of salience)
        epsMat = eye(size(TF,1))* 1e-12 ; 
        E = ones(numFrame,1) ; 
        O = zeros(numFrame,1) ;
        U = eye(numFrame) ;     
        
%         nRows = numFrame*(numFrame-1)/2  ; 
%         A = zeros( nRows,numFrame) ; 
%         cnt = 1 ; 
%         for p = 1 : numFrame
%             for q = p + 1 : numFrame
%                 A(cnt,p) = 1 ; 
%                 A(cnt,q) = -1 ; 
%                 cnt = cnt + 1 ; 
%             end 
%         end
        
        nRows = numFrame ; 
        A = zeros( nRows,numFrame) ; 
        for p = 1 : numFrame-1
            A(p,p) = 1 ; 
            A(p,p+1) = -1 ; 
        end

        M = w'*TX*TL ; 
        N = w'*TX*TL' ; 
        H11 = M'*M + gamma * TF + yita*U + epsMat ; 
        H12 = zeros(numFrame) ; 
        H21 = zeros(numFrame) ; 
        H22 = augRate * N'*N + gamma * TF2 + yita*U + epsMat ; 
        H = 2*[H11, H12; H21, H22] ; 
        f = [- 2 * TY * M' ; -2 * augRate * TY0 * N'] ; 
        Aeq = [E'*TL,O';O',E'*TL'; A*TL, A*TL'] ; 
        beq = [1;1; zeros(nRows,1)] ; 
        lb = [O;O] ; 
        TB0 = [beta{i}; betaS{i}]; 
        
        if ~issymmetric(H)
           H = (H'+H)/2 ; 
        end
        TBB = quadprog(H,f,[],[],Aeq,beq, lb,[],TB0,options);
        
        beta{i} = TBB(1:numFrame) ; 
        alpha{i} = TL*beta{i} ; 
        

        betaS{i} = TBB(numFrame+1:end) ; 
        alphaS{i} = TL'*betaS{i} ;
        % ===============================================
    end

    temScore = funScore(w,alpha,alphaS,beta,betaS,LapMat,Y,Y0,X,FMat,...
                            lambda,gamma,rho,isWeighted,augRate) ; 
    score = [score;temScore] ; 

    fprintf('== Iter %d , Score: %.3f -> %.3f ...\n', TT,score(end-1),score(end)) ; 

    WRes = [WRes,{w}] ; 
    alphaRes = [alphaRes;{alpha}] ; 
    alphaResS = [alphaResS;{alphaS}] ; 

    if abs(score(end) - score(end-1)) < thresh
        break; 
    end
end


function [score] = funScore(w,alpha,alphaS,beta,betaS,LapMat,Y,Y0,X,FMat,...
                            lambda,gamma,rho,isWeighted,augRate)
%% compute the objective function score 

OX = cell(size(alpha)) ; % average instance
OXS = cell(size(alphaS)) ; % average instance
for i = 1 : length(OX)
    OX{i} = X{i}*alpha{i} ; 
    OXS{i} = X{i}*alphaS{i} ; 
end
[OX,numVec] = cellToMat(OX) ; 
[OXS,numVec] = cellToMat(OXS) ;
OY = Y(:) ; 
OYS = Y0(:); 

sampW = ones(length(OY),1) ;
if isWeighted
    numNeu = sum(OY == 0); 
    numPos = length(OY) - numNeu ; 
%     sampW(OY == 0) = sqrt(1 / numNeu) ; 
%     sampW(OY ~= 0) = sqrt(1 / numPos) ; 
    sampW(OY == 0) = 1 ; 
    sampW(OY ~= 0) = sqrt(numNeu / numPos) ; 
end

OD = diag(sampW.*sampW) ;
A = (OY - OX' * w)' * OD * (OY - OX' * w) ; 
A2 = (OYS - OXS' * w)' * OD * (OYS - OXS' * w) ;
B =  w' * w ; 
C = w' * FMat * w ; 


DA = 0 ; 
DB = 0 ; 
for i =1 : length(LapMat)
    TF = LapMat{i} ;
    
    TA = alpha{i} ; 
    TAS = alphaS{i} ;
    DA = DA + TA' * TF * TA  + TAS' * TF * TAS; 
    
    TB = beta{i} ; 
    TBS = betaS{i} ;
    DB = DB + TB' * TF * TB  + TBS' * TF * TBS; 
end

score = A + augRate * A2 + lambda * B  + rho * C + gamma * DA ; 



