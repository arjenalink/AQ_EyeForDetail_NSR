%% This script imports the RAW data (not publically available) and saves it in a more compact manner
% in 'SelFts_Sj_And_Beh_Data.mat'. Furthermore, it creates the ImportAndSaveAllData_Master_DATA_WithRTsNOL.mat
% file used by the main data analysis script
projDir='C:\Users\alink\Google Drive\diagnosticFeatures\experiment_v2\';
SjFolders=dir([projDir '\']);
%% global variables
nAQquestions=50;
nSubScales=5;
nIms=10;
FtPerIm=1000;
nTrialsPerRun=500;
nFtPresented=90;
%% Select the set of subjects for which both AQ and Visual data is present 
[NUM,TXT,RAW]=xlsread('AnswersToAllQuestionsInCorrSjOrder.xlsx');
AllSjNames=cell(size(RAW,2),1);
for i=1:numel(AllSjNames)
    if isnumeric(RAW{1,i})
        AllSjNames{i}=num2str(RAW{1,i});
    else
        AllSjNames{i}=RAW{1,i};
    end
end
FldNames={};
for f=1:numel(SjFolders)
    FldNames{f}=SjFolders(f).name;
end
SjWithDataRef=ismember(AllSjNames,FldNames);
SjWithDataRef(strcmp(AllSjNames,'lm020898')+strcmp(AllSjNames,'71467ALD')==1)=0; % two subjects with corrupted data
SjWithDataNames=AllSjNames(SjWithDataRef==1);
nSjs=numel(SjWithDataNames);

%% Save all AQ data per sj
All50AQanswersPerValidSj=NUM(2:nAQquestions+1,SjWithDataRef); %first line contains sj names
[NUM,TXT,RAW]=xlsread('QuestionWeights_PerSubscore.xlsx');
Q_AQ_weights=sum(NUM,2);
QsSubscales=zeros(nAQquestions,1);
for i=1:nAQquestions
    QsSubscales(i)=find(NUM(i,1:nSubScales)~=0);
end
WeigthedAnswers=All50AQanswersPerValidSj.*repmat(Q_AQ_weights,1,nSjs);
TraitPresent=WeigthedAnswers==1|WeigthedAnswers==2|WeigthedAnswers==-3|WeigthedAnswers==-4;
AqScoresFromRawData=sum(TraitPresent,1);
temp=TraitPresent.*repmat(QsSubscales,1,nSjs);
subAQscores=zeros(nSjs,nSubScales);
for sc=1:nSubScales
    subAQscores(:,sc)=sum(temp==sc,1);
end
[NUM,TXT,RAW]=xlsread('PercHavingTrait_AS_Ctr_Students.xlsx');
TraitPrevalence_AQ_CONTR_STUDENTS=NUM(1:nAQquestions,2:4);
QuestionClDiag=NUM(1:nAQquestions,2)./NUM(1:nAQquestions,4);

%% Read in all behvioural data in same order as the AQ sj order
FtSel=zeros(0,nFtPresented);
SjRef=[];
BehDat=[];
for sj=1:nSjs
    sj
    tic
    BehDatFiles=dir([projDir SjWithDataNames{sj} '\*.txt']);
    FtInfFiles=dir([projDir SjWithDataNames{sj} '\*.mat']);
    if numel(BehDatFiles)==numel(FtInfFiles)        
        for run=1:numel(BehDatFiles)
            if run==1&sj==1
                temp=importdata([projDir SjWithDataNames{sj} '\' BehDatFiles(run).name]);
                BehDat=temp.data;
                FtInf=load([projDir SjWithDataNames{sj} '\' FtInfFiles(run).name]);
                for i=1:nTrialsPerRun
                    FtSel(end+1,:)=find(FtInf.params.SelFtsPerTrial(i,:)==1)';
                end
                SjRef=[SjRef;ones(nTrialsPerRun,1)*sj];
            else
                FtInf=load([projDir SjWithDataNames{sj} '\' FtInfFiles(run).name]);
                for i=1:nTrialsPerRun
                    FtSel(end+1,:)=find(FtInf.params.SelFtsPerTrial(i,:)==1)';
                end
                temp=importdata([projDir SjWithDataNames{sj} '\' BehDatFiles(run).name]);
                BehDat=[ BehDat; temp.data];
                SjRef=[SjRef;ones(nTrialsPerRun,1)*sj];
            end            
            toc
        end
    end
end
BehDatColumNames={'trialnumber',	'ImageNr',	'ImageCategory',	'trialOnset',	'rt',	'ImageCategory_Response_3isDontKnow',	'runNo'};	
save('SelFts_Sj_And_Beh_Data.mat','BehDat','SjRef','FtSel','BehDatColumNames');
%% Compute prop accuracy per feature
FtDiagnRatingPerSj=zeros(nSjs,FtPerIm,nIms);
nObsPerFtPerSj=zeros(nSjs,FtPerIm,nIms);
FtRTPerSj=zeros(nSjs,FtPerIm,nIms);
FtRTPerSjNOL=zeros(nSjs,FtPerIm,nIms);
RTsPerSjIm=zeros(nSjs,nIms);
PerfPerSjIm=zeros(nSjs,nIms);
RTsPerSjImNOL=zeros(nSjs,nIms);
for sj=1:nSjs
    sj
 
     thisSj=SjRef==sj;
    for im=1:nIms
          
        CorrRef=thisSj&BehDat(:,2)==im&BehDat(:,3)==BehDat(:,6);
        PerfPerSjIm(sj,im)=mean(CorrRef)/mean(thisSj&BehDat(:,2)==im);
        
        
        RTs=BehDat(CorrRef,5);
        RTSzscored=zscore( RTs);
        RTsPerSjIm(sj,im)=mean(RTs);
        RTsPerSjImNOL(sj,im)=mean(RTs(abs(RTSzscored)<2));
        for ft=1:FtPerIm
            FtPresRef=thisSj&sum(FtSel==ft,2)==1;
            CorrRef=BehDat(:,2)==im&BehDat(:,3)==BehDat(:,6)&FtPresRef;
            AccRate=mean(CorrRef)/mean(thisSj&BehDat(:,2)==im&FtPresRef);
            RTs=BehDat(CorrRef,5);
            RTSzscored=zscore( RTs);
            FtRTPerSj(sj,ft,im)=mean(RTs);
            FtRTPerSjNOL(sj,ft,im)=mean(RTs(abs(RTSzscored)<2));
            FtDiagnRatingPerSj(sj,ft,im)=AccRate;
            nObsPerFtPerSj(sj,ft,im)=sum(FtPresRef&BehDat(:,2)==im);
        end
    end  % acsassccc SelFtMem=[]; SelFtMem=[]; SelFtMem=[]; SelFtMem=[]; SelFtMem=[];
    toc
end
%% save all relevant data in one .mat file
save('ImportAndSaveAllData_Master_DATA_WithRTsNOL.mat','FtDiagnRatingPerSj','nObsPerFtPerSj','AqScoresFromRawData','QuestionClDiag','TraitPresent','WeigthedAnswers','subAQscores','Q_AQ_weights','QsSubscales','FtInf','FtRTPerSj','FtRTPerSjNOL','RTsPerSjIm','RTsPerSjImNOL','PerfPerSjIm','-v7.3');
