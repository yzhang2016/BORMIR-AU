function [trX,trY,trY0,tsX,tsY,trGTY] = prepareBP4D_AU(seqs,DatInd,varargin)
%% preparing data for BP4D AU 

featStr = 'LMark' ; 
for i =1 :2: length(varargin)
    if strcmp(varargin{i},'feature')
        featStr = varargin{i+1} ; 
    else
        error('Not such a key...') ; 
    end
end

Feat = [] ; 
AUINT = [] ; 
SUBInd = [] ; 

for i = 1:length(seqs)
    tem = seqs(i) ; 
    if strcmp(featStr,'NormLBP')
        Feat = [Feat;{tem.LBPFeat(:,1:end-1)}] ;  % LBP feature (97% energy)
    elseif strcmp(featStr,'PCANormLMark')
        Feat = [Feat;{tem.LMark_PCA(:,1:end-1)}] ; % PCA landmark fetaure (97% energy)
    elseif strcmp(featStr,'NormLMark')
        Feat = [Feat;{tem.LMark_Norm(:,1:end-1)}] ; % nomarlized landmark fetaure 
    else
        error('Wrong: not such feature') ; 
    end

    AUINT = [AUINT; {tem.AUINT}] ; 
    SUBInd = [SUBInd;tem.SUBInd] ; 
end

trainInd = DatInd.trainInd  ; 
testInd = DatInd.testInd ; 

trainFeat = [] ; 
trainINT = [] ; 
testFeat = [] ; 
testINT = [] ; 

numTrSeq  = length(trainInd) ; 
numTsSeq = length(testInd) ; 

%% testing 
for  i = 1 : numTsSeq
    temInd = testInd(i) ; 
    temFeat = Feat{temInd}'; 
    tem = AUINT{temInd}(:,1)' ; 
    temINT = AUINT{temInd}(:,2)';

    if tem(1) > tem(end)
        temFeat = temFeat(:,end:-1:1) ; 
        temINT = temINT(end:-1:1) ; 
    end
    
    testFeat = [testFeat,temFeat] ; 
    testINT = [testINT,temINT] ; 
end

tsX = testFeat; 
tsY = testINT ; 

%% training 
trY = [] ; 
trY0 = [] ; 
trGTY = [] ; 
for  i = 1 : numTrSeq
    temInd = trainInd(i) ; 
    temFeat = Feat{temInd}'; 
    temINT = AUINT{temInd}(:,2)' ; 
    
    trainFeat = [trainFeat,{temFeat}] ; 
    trGTY = [trGTY ; temINT']; 
    trY = [trY,temINT(end)] ; 
    trY0 = [trY0,temINT(1)] ; 
end
trX = trainFeat ; 

%% to keep consistent 
numSel = floor(sum(DatInd.selInd)/2) ; 

if numSel > length(trX)
    numSel = length(trX) ; 
end

trX = trX(1:numSel); 
trY = trY(1:numSel) ; 
trY0 = trY0(1:numSel) ; 




