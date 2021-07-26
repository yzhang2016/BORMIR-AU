clear; clc ; 
% delete(gcp('nocreate')) ; 
warning off ; 
addpath('../data') ; 
addpath('./utility') ; 
addpath('../') ; 


%% Configuration for BP4D
dataName = 'BP4D'; % BP4D
dataIndName = 'BP4D' ; % BP4D
featType = 'NormLMark' ; 
AUInd = [6,10,12,14,17] ; % BP4D 
rateRange = 1 ;

for ind = 1: length(AUInd)
%% data process
dataPath = sprintf('../data/%s/AU/AUData_lmark_AU%d.mat',dataName,AUInd(ind)) ; 
src = load(dataPath) ; 
seqs = src.seqs; 
cvPath = sprintf('./%s_AU_5fds_protol/AU%d',dataIndName,AUInd(ind));
dstPath = sprintf('../Result_AU_protol/%s/BOMIR/AU%d',dataIndName,AUInd(ind)) ; 
if ~exist(dstPath,'dir')
    mkdir(dstPath) ; 
end

option.isRmvOrdinal = 0 ; 
option.maxIter = 20 ; 
option.thresh = 1e-3 ;
option.isWeighted = 0 ;
option.isEqualFrameW = 1; 


%COMMON Parameters 
% As the data distribution of each AU differs, these parameters 
% need to be tuned.
augRate = 0.5;  % augmented
lambda = 5;    % L2 on W
rho =  0.3;    % intensity smoothness
gamma =  1;  % ordinal smoothness 
yita = 1;     % a1+b1 =a2 + b2 


%% leave one subject out 
numTimes = 5; 

for JJ = 1 : length(rateRange)
annoRate = rateRange(JJ);

TT_tsRES = [] ; 

cvname = sprintf('%s/AnnoRate_%.2f.mat',cvPath,annoRate) ; 
cvDat = load(cvname); 
cvDat = cvDat.oneRate ; 

if annoRate <= 0.1 || annoRate == 1
    TemNumTimes = 1 ; 
else
    TemNumTimes = numTimes ; 
end

for TT = 1 : TemNumTimes 

TTCvDat =  cvDat{TT} ; 

numFds = length(TTCvDat) ; 
tsRES =  zeros(numFds,4) ; 

%par
for i =1 : numFds
    
    SUBDatInd = TTCvDat(i) ; 
    
    % data prepare
    [trX,trY,trY0,tsX,tsY,trGTY] = prepareBP4D_AU(seqs,SUBDatInd,'feature',featType) ;   
     
    tic ; 
    [w,alpha,beta,score,WRes,alphaRes,alphaResS] =...
            BOMIR_IV_B_ext(trY,trY0,trX,lambda,gamma,rho,yita,augRate,option);
     
    endtime = toc ;      

    tic 
    [predVal] = SequenceTest(w,tsX) ; 
    endtime2 = toc ; 
    
    [PCC,UICC,UMAE,UMSE] = OSWMeasure(predVal,tsY)  ; 
    
    tsRES(i,:) = [PCC,UICC,UMAE,UMSE] ; 
    
    
    %%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%
    figure ; 
    pos = alphaRes{end} ; 
    neg = alphaResS{end} ; 
    for MM = 1 : 16
        temA = pos{MM}' ; 
        temA2 = neg{MM}' ; 
        subplot(4,4,MM) ;

        plot([temA',temA2'],'LineWidth',2) ; 
        xlim([1,length(temA)]) ; 
        ylim([0,0.5]) ; 
        tl = ['Intensity [',num2str(trY(MM)),']'] ;
        title(tl); 
    end 
    
    fprintf('Rate = %.2f, Time = %d, SUB = %d...\n',annoRate,TT,i); 
end

avgTSRES = mean(tsRES,1) ; 
TT_tsRES = [TT_tsRES;avgTSRES] ; 

end

avgTT_tsRES = mean(TT_tsRES,1) ; 
stdTT_tsRES = std(TT_tsRES,1) ; 

svname = sprintf('%s/rate_%.2f.mat',dstPath,annoRate) ; 
save(svname,'avgTT_tsRES','stdTT_tsRES','TT_tsRES') ; 
end

end