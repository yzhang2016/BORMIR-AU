function [outTrFeat,outY,outTsFeat,outTsLb] = dataPreprocessCKBP4D(trFeat,trLb,tsFeat,tsLb,numDim)
%% processing 
% 1. use the first and the last few frames to synthesize a sequence 
% 

negFrameNum = 0; 16;

% Use both onset and downset parts
outTrFeat = [] ; 
outY = [] ; 

temFeat = cell2mat(trFeat) ; 
temFeat = [temFeat,tsFeat] ; 
temFeat = temFeat' ; 

coeff = pca(temFeat) ; 
numComp = numDim ; 
coeff = coeff(:,1:numComp) ; 

for i =1 : length(trFeat)
    temM = trLb{i} ; 
    
    if temM(end,2) >= temM(end-1,2)
        outY = [outY,temM(end,2)] ; 
        outTrFeat = [outTrFeat, {coeff'*trFeat{i}}] ; 
        
         % to synthesize sequence with the intensity 0
         startNum = 1 ; % only use the first frames and the last frames
         repTimes = floor(negFrameNum / (startNum)) ; 
          
        sFeat = trFeat{i}(:,1:startNum) ; 
        eFeat = [] ; %trFeat{i}(:,end:end-startNum+1); 
        cFeat = [sFeat,eFeat] ; 
        cFeat = repmat(cFeat,1,repTimes) ; 
        cY = 0 ; 
        outY = [outY,cY] ; 
        outTrFeat = [outTrFeat,{coeff'*cFeat}] ; 

    else
        outY = [outY,temM(end-1,2)] ;
        outTrFeat = [outTrFeat,{coeff'* (trFeat{i}(:,1:temM(end-1,1)))}];
        
        outY = [outY,temM(end-1,2)] ;
        outTrFeat = [outTrFeat,{coeff'*(trFeat{i}(:,end:-1:temM(end-1,1)))}];
         % to synthesize sequence with the intensity 0
         startNum = 1 ; % only use the first frames and the last frames 
         repTimes = ceil(negFrameNum / (startNum * 2 )) ; 
         
        sFeat = trFeat{i}(:,1:startNum) ; 
        eFeat = sFeat ; %trFeat{i}(:,end:-1:end-startNum+1); 
        cFeat = [sFeat,eFeat] ; 
        cFeat = repmat(cFeat,1,repTimes) ; 
        cY = 0 ; 
        outY = [outY,cY] ; 
        outTrFeat = [outTrFeat,{coeff' *cFeat}] ; 
    end

end

abandon = [] ; 
for i = 1: length(outTrFeat)
    if isempty(outTrFeat{i})
        abandon = [abandon,i] ; 
    end
end
outTrFeat(abandon)= [] ; 
outY(abandon) = [] ; 

if isempty(tsFeat) || isempty(tsLb)
    outTsFeat = [] ; 
    outTsLb = [] ; 
else
    outTsFeat = coeff' *tsFeat ;
    outTsLb = tsLb ; 
end
