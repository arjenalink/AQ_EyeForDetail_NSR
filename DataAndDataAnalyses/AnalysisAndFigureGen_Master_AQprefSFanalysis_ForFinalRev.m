%% This script performs all data analyses and generates all results reported 
% in our Nature Scientific Reports paper "Clinically relevant autistic
% traits predict greater reliance on detail for image recognition", 2020 
clear all;close all
%% load .mat files containing all behavioural data
load ImportAndSaveAllData_Master_DATA_WithRTsNOL.mat
load SelFts_Sj_And_Beh_Data.mat
%% Global variables
nSjs=size(FtDiagnRatingPerSj,1);
nAQquestions=numel(Q_AQ_weights);
nSubScales=5;
IMSIze=250;
SelExemplImage=7;
Displ=zeros(250);
nFtsPerIm=1000;
NrIms=10;
nFtPropConsidered=4;
NrBins=5;
CatRef=1:5;
DogRef=6:10;
ImSizeInVDA=22.5;
%% load feature parameters for all 10 images and save them in one FtParamMatrix
cyclesPerVOF=squeeze(FtInf.params.ImGabParams(:,2,:));
DistFromCentre=squeeze(sqrt((abs((FtInf.params.ImGabParams(:,5,:)-.5)).^2)+(abs((FtInf.params.ImGabParams(:,6,:)-.5)).^2)));
GabOris=squeeze(FtInf.params.ImGabParams(:,3,:));
load MinDistFromEye.mat;
SF_DFE_DFC_Ori_FtsParams=zeros(nFtPropConsidered,NrIms*nFtsPerIm);
SF_DFE_DFC_Ori_FtsParams(1,:,:)=cyclesPerVOF(:);
SF_DFE_DFC_Ori_FtsParams(2,:,:)=MinDistFromEye(:);
SF_DFE_DFC_Ori_FtsParams(3,:,:)=DistFromCentre(:);
SF_DFE_DFC_Ori_FtsParams(4,:,:)=GabOris(:);
ImInds=repmat(1:10,1000,1);ImInds=ImInds(:);
CatDogInds=(ImInds<6)+1;
%% remove singel NaN values and z-transform FDi values and vectorize across images
FtDiagnRatingPerSj(isnan(FtDiagnRatingPerSj(:)))=.5;
FtDiagnRatingPerSjNorm=zscore(FtDiagnRatingPerSj,[],2);
FtDiagnRatingPerSjNormVect=FtDiagnRatingPerSjNorm(:,:);

%% prepare reaction time data
FtRTPerSj(isnan(FtRTPerSj))=mean(mean(FtRTPerSj(isnan(FtRTPerSj)==0)));
ResponseSpeedNorm=-zscore(FtRTPerSj,[],2);
ResponseSpeedNormVect=ResponseSpeedNorm(:,:);

%% Define high and low AQ groups
medAQ=median(AqScoresFromRawData);
LowAQsjs=find(AqScoresFromRawData<=medAQ);
HighAQsjs=find(AqScoresFromRawData>medAQ);

%%  behavioural performance results
%nObsWeights=nObsPerFtPerSj(:)/mean(nObsPerFtPerSj(:));
GrandMean=mean(PerfPerSjIm(:))
PerfPerSj=mean(PerfPerSjIm,2);
PerfSTD=std(PerfPerSj)
PerfPerSj_CatDog=[mean(PerfPerSjIm(:,CatRef),2),mean(PerfPerSjIm(:,DogRef),2)] ;

AQgroup=repmat([ones(numel(LowAQsjs),1);ones(numel(HighAQsjs),1)*2],1,2);
CatDog=[ones(nSjs,1),ones(nSjs,1)*2];
Cat_Dog_Perf=mean(PerfPerSj_CatDog)
[h p ci stats]=ttest(PerfPerSj_CatDog(:,1)-PerfPerSj_CatDog(:,2));
t= stats.tstat
anovan(PerfPerSj_CatDog(:),[AQgroup(:),CatDog(:)],'full')
LAQperf=mean(mean(mean(FtDiagnRatingPerSj(LowAQsjs,:,:)),2),3)
HAQperf=mean(mean(mean(FtDiagnRatingPerSj(HighAQsjs,:,:)),2),3)

ProportionUnsurePerSj=zeros(nSjs,1);
ProportionIncorrectPerSj=zeros(nSjs,1);
for sj=1:nSjs
    ProportionUnsurePerSj(sj)=mean(BehDat(SjRef==sj,6)==3);
    ProportionIncorrectPerSj(sj)=mean(BehDat(SjRef==sj,6)<3&BehDat(SjRef==sj,3)~=BehDat(SjRef==sj,6));
end
ProportionUnsurePerSj_GrandMean=mean(ProportionUnsurePerSj)
ProportionUnsurePerSj_SD=std(ProportionUnsurePerSj)

ProportionIncorrectPerSj_GrandMean=mean(ProportionIncorrectPerSj)
ProportionIncorrectPerSj_SD=std(ProportionIncorrectPerSj)
%%  Behavioural RT results
GrandMean=nanmean(RTsPerSjIm(:))
RTPerSj=nanmean(RTsPerSjIm(:,:),2);
PerfSTD=std(RTPerSj)
RTPerSj_CatDog=[mean(mean(FtRTPerSj(:,:,CatRef),3),2),mean(mean(FtRTPerSj(:,:,DogRef),3),2)] ;
AQgroup=repmat([ones(numel(LowAQsjs),1);ones(numel(HighAQsjs),1)*2],1,2);
CatDog=[ones(nSjs,1),ones(nSjs,1)*2];
Cat_Dog_RT=mean(RTPerSj_CatDog)
[h p ci stats]=ttest(RTPerSj_CatDog(:,1)-RTPerSj_CatDog(:,2));
t= stats.tstat
anovan(PerfPerSj_CatDog(:),[AQgroup(:),CatDog(:)],'full')
LAQRT=mean(RTPerSj(LowAQsjs))
HAQRT=mean(RTPerSj(HighAQsjs))
[h p ci stats]=ttest2(RTPerSj(LowAQsjs),RTPerSj(HighAQsjs))

%% FIGURE 2
% RFX anova test for AQ interaction with SF, DFE,DFC and Ori effects
AQref=zeros(nSjs,1);AQref(LowAQsjs)=1;AQref(HighAQsjs)=2;
MeanMinMaxVals=zeros(nFtPropConsidered,5,3);
AllFDiVals=zeros(nFtPropConsidered,nSjs*NrBins,1);
AllBinAQpreds=zeros(nFtPropConsidered,nSjs*NrBins,2);
figure(99);
DimNames={'SF','DFE','DFC'};
for dim=1:numel(DimNames)
    FDiVals=zeros(nSjs*NrBins,1);
    BinAQpreds=zeros(nSjs*NrBins,2);
    binSEMs=zeros(NrBins,2);
    binMeans=zeros(NrBins,2);
    [vals sortedFtInds]=sort(SF_DFE_DFC_Ori_FtsParams(dim,:));
    BinInds=ceil(([1:nFtsPerIm*NrIms]*NrBins)/(nFtsPerIm*NrIms));
    count=1;
    for sj=1:nSjs
        for b=1:NrBins
            BinFts=sortedFtInds(BinInds==b);
            FDiVals(count)=mean(FtDiagnRatingPerSjNormVect(sj,BinFts));
            BinAQpreds(count,:)=[b, AQref(sj)];
            if dim==1
                MeanMinMaxVals(dim,b,:) = [mean(vals(BinInds==b)) ,min(vals(BinInds==b)),max(vals(BinInds==b))]./ImSizeInVDA;
            else
                MeanMinMaxVals(dim,b,:) = [mean(vals(BinInds==b)) ,min(vals(BinInds==b)),max(vals(BinInds==b))]*ImSizeInVDA;
            end
            if sj==1
                Displ(:)=squeeze(mean(FtInf.params.gabvects(:,BinFts(ismember(BinFts,find(ImInds==SelExemplImage)))),2));
                figure(88); imagesc(Displ);colormap('gray');axis square;set(gca,'xtick',[]);set(gca,'ytick',[]); set(gca,'LooseInset',get(gca,'TightInset'));saveas(gcf,[DimNames{dim} '_Bin' num2str(b) 'MeanMinMax' num2str(MeanMinMaxVals(dim,b,:)) '.png']);
            end
            count=count+1;
        end
    end
    for b=1:NrBins
        binSEMs(b,1)= std(FDiVals(BinAQpreds(:,1)==b&BinAQpreds(:,2)==1))/sqrt(sum(AQref==1));
        binSEMs(b,2)= std(FDiVals(BinAQpreds(:,1)==b&BinAQpreds(:,2)==2))/sqrt(sum(AQref==2));
        binMeans(b,1)= mean(FDiVals(BinAQpreds(:,1)==b&BinAQpreds(:,2)==1));
        binMeans(b,2)= mean(FDiVals(BinAQpreds(:,1)==b&BinAQpreds(:,2)==2));
    end
    figure(99);
    subplot(1,3,dim);hold on;
    errorbar(binMeans(:,1),binSEMs(:,1),'linewidth',3,'color',[.1 .1 .1; ]) ; xlim([.5 5.5]); ylim([-.02 .038]); ax=gca; set(ax,'Xticklabels',[]);set(ax,'yticklabels',[]);set(ax,'LineWidth',1.5);
    errorbar(binMeans(:,2),binSEMs(:,2),'linewidth',3,'color',[.8 .1 .1; ]) ; xlim([.5 5.5]);ylim([-.02 .038]); ax=gca; set(ax,'Xticklabels',[]);set(ax,'yticklabels',[]);set(ax,'LineWidth',1.5);
    AllFDiVals(dim,:,:)=FDiVals;
    AllBinAQpreds(dim,:,:) =BinAQpreds;
end
% plot selected anova results, 1:3={'SF','DFE','DFC'}
selDim=3;
figure(222);
anovan(AllFDiVals(selDim,:,:)',squeeze(AllBinAQpreds(selDim,:,:)),'full');

%% FIGURE 2
% assess how well individual traits predict SF preference interaction effect size
%(Make sure this is ran AFTER the above section and not the RT section
FDiPerBinAndSj=zeros(52,NrBins);
for sj=1:52
    SjInds=[1:NrBins]+(NrBins*(sj-1))
    FDiPerBinAndSj(sj,:)=AllFDiVals(1,SjInds);
end
[h p ci stats]=ttest2(FDiPerBinAndSj(HighAQsjs,5),FDiPerBinAndSj(LowAQsjs,5));

B=zeros(nSjs,1);
for sj=1:nSjs
    temp=regress(FDiPerBinAndSj(sj,:)',[1 2 3 4 5 ; 1 1 1 1 1]');
    B(sj)=temp(1,1);
end
SFmod=B;

SFprefModPerTrait=zeros(50,1);
for q=1:nAQquestions
    SFprefModPerTrait(q)= mean(SFmod(TraitPresent(q,:)==1))-mean(SFmod(TraitPresent(q,:)==0));
end

ScInds=repmat(1:5,10,1);
[p,t,stats]=anova1(SFprefModPerTrait,QsSubscales);
[c,m,h,nms] = multcompare(stats,'cType','scheffe');
[vals, inds]=sort(SFprefModPerTrait,'descend');
% five  posthoc t-tests if a trait type gives rise to greater or smaller contribution to H-SF preference
PerSCdata=[SFprefModPerTrait(QsSubscales==1),SFprefModPerTrait(QsSubscales==2),SFprefModPerTrait(QsSubscales==3),SFprefModPerTrait(QsSubscales==4),SFprefModPerTrait(QsSubscales==5)];

for sc=1:5
    sc
    OtherVals=PerSCdata(:,setxor(sc,1:5));
    [h p ci stats]=ttest2(PerSCdata(:,sc),OtherVals(:));
    stats.tstat
    p=p*5 %bonf correction
end
OtherVals=PerSCdata(:,[2 4 5]);
[h p ci stats]=ttest2(PerSCdata(:,1),OtherVals(:));
OtherVals=PerSCdata(:,[2 4 5])
[h p ci stats]=ttest2(PerSCdata(:,3),OtherVals(:));

% figure(78);notBoxPlot(SFmodPerTraitType,1:5,'interval' ,'tInterval');
figure(77);
barwitherr(std(PerSCdata)/sqrt(10),mean(PerSCdata))
ylim([-.001 .005]);
set(gca,'YTick',-.001:.001:.004)
set(gca,'Xticklabel',[])
set(gca,'Xtick',[])


[B,STATS] =robustfit(log(QuestionClDiag),SFprefModPerTrait);
STATS.t
STATS.p
figure(55); hold on
[vals inds]=sort(log(QuestionClDiag));
plot(vals,B(1)+[B(2)*log(QuestionClDiag(inds))],'r','LineWidth',3);
scatter(log(QuestionClDiag),SFprefModPerTrait,50,'k','filled');
xlabel('item clinical diagnosticity [log odds ratio)]');ylabel('item-realted high spatial frequency preference increase [a.u.]');
[rho p]=corr(log(QuestionClDiag),SFprefModPerTrait);


%% FIGURE 1c
% check test retest reliability of FDis by assesing the sign. of the split half correlation
% set ComputeNow=1 to recompute (takes a while), and 0 to load precomputed  data
ComputeNow=0;
nCorsPerAv=100;
if ComputeNow==1
    nPerms=10000;
else
    nPerms=1;
end

shCorrVals=zeros(nCorsPerAv,1);
rng('default');
rng(42);
tic
for i=1:nCorsPerAv
    SelDat=FtDiagnRatingPerSjNormVect(randperm(nSjs),:);
    MeanOddVect=mean(SelDat(1:2:end,:));
    MeanEvenVect=mean( SelDat(2:2:end,:));
    shCorrVals(i)=corr(MeanOddVect',MeanEvenVect');
end
toc
obsSHcorr=mean(shCorrVals);
shCorrValsPermed=shCorrVals;
meanPermRhos=zeros(nPerms,1);
for perm=1:nPerms
    perm
    tic
    SelDat=FtDiagnRatingPerSjNorm;
    for sj=1:nSjs
        for im=1:10
            SelDat(sj,:,im)=FtDiagnRatingPerSjNorm(sj,randperm(1000),im);
        end
    end
    for i=1:nCorsPerAv
        SelDat=SelDat(randperm(48),:);
        MeanOddVect=mean( SelDat(1:2:end,:));
        MeanEvenVect=mean( SelDat(2:2:end,:));
        shCorrValsPermed(i)=corr(MeanOddVect',MeanEvenVect');
    end
    meanPermRhos(perm)=mean(shCorrValsPermed);
    toc
end

if ComputeNow==0
    load ObservedAndPermedAvSHcorrVals;
end
[n, x]=hist(meanPermRhos,80);
hist(meanPermRhos,80);
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [0.5 0.5 0.5];
line([mean( shCorrVals) obsSHcorr], [0 max(n)],'Color',[.8 .1 .1],'LineWidth',3);
xlim([-.03 .09])
ax=gca;
ax.YTick= 0:100:500;
ax.LineWidth=2;
set(gca, 'TickDir', 'out');

%% RFX anova test for AQ interaction with SF, DFE,DFC and Ori effects, REACTION TIMES
AQref=zeros(nSjs,1);AQref(LowAQsjs)=1;AQref(HighAQsjs)=2;
MeanMinMaxVals=zeros(nFtPropConsidered,5,3);
AllFDiVals=zeros(nFtPropConsidered,nSjs*NrBins,1);
AllBinAQpreds=zeros(nFtPropConsidered,nSjs*NrBins,2);
figure(199);
DimNames={'SF','DFE','DFC'};
for dim=1:numel(DimNames)
    FDiVals=zeros(nSjs*NrBins,1);
    BinAQpreds=zeros(nSjs*NrBins,2);
    binSEMs=zeros(NrBins,2);
    binMeans=zeros(NrBins,2);
    [vals sortedFtInds]=sort(SF_DFE_DFC_Ori_FtsParams(dim,:));
    BinInds=ceil(([1:nFtsPerIm*NrIms]*NrBins)/(nFtsPerIm*NrIms));
    count=1;
    for sj=1:nSjs
        for b=1:NrBins
            BinFts=sortedFtInds(BinInds==b);
            FDiVals(count)=mean(ResponseSpeedNormVect(sj,BinFts));
            BinAQpreds(count,:)=[b, AQref(sj)];
            MeanMinMaxVals(dim,b,:) = [mean(vals(BinInds==b)) ,min(vals(BinInds==b)),max(vals(BinInds==b))];
            count=count+1;
        end
    end
    for b=1:NrBins
        binSEMs(b,1)= std(FDiVals(BinAQpreds(:,1)==b&BinAQpreds(:,2)==1))/sqrt(sum(AQref==1));
        binSEMs(b,2)= std(FDiVals(BinAQpreds(:,1)==b&BinAQpreds(:,2)==2))/sqrt(sum(AQref==2));
        binMeans(b,1)= mean(FDiVals(BinAQpreds(:,1)==b&BinAQpreds(:,2)==1));
        binMeans(b,2)= mean(FDiVals(BinAQpreds(:,1)==b&BinAQpreds(:,2)==2));
    end
    figure(199);
    subplot(1,3,dim);hold on;
    errorbar(binMeans(:,1),binSEMs(:,1),'linewidth',3,'color',[.1 .1 .1; ]) ;  ax=gca; set(ax,'Xticklabels',[]);set(ax,'LineWidth',1.5);
    errorbar(binMeans(:,2),binSEMs(:,2),'linewidth',3,'color',[.8 .1 .1; ]) ;  ax=gca; set(ax,'Xticklabels',[]);set(ax,'LineWidth',1.5);
    AllFDiVals(dim,:,:)=FDiVals;
    AllBinAQpreds(dim,:,:) =BinAQpreds;
end
% plot selected anova results
selDim=3;
figure(2222);
anovan(AllFDiVals(selDim,:,:)',squeeze(AllBinAQpreds(selDim,:,:)),'full');
%[c,m,h,nms] = multcompare(stats,'CType','lsd')
% checking if AQxSF interaction depends on reaction times
[rho p]=corr(SFmod,nanmean(FtRTPerSj(:,:),2))

rts=nanmean(FtRTPerSj(:,:),2);
[rho p]=corr(SFmod,rts);
scatter(SFmod,rts);
[B,STATS] =robustfit(SFmod,rts);
[vals inds]=sort(SFmod);
figure(1000); hold on
plot(vals,B(1)+[B(2)*SFmod(inds)],'r','LineWidth',3);
scatter(SFmod,rts,50,'k','filled')
xlabel('spatial frequency preference modulation');ylabel('reaction times');

%% check test retest reliability of RT-FDis by assesing the sign. of the split half correlation
% set ComputeNow=1 to recompute (takes a while), and 0 to load precomputed  data
ComputeNow=0;
nCorsPerAv=100;
if ComputeNow==1
    nPerms=10000;
else
    nPerms=1;
end

nCorsPerAv=100;
shCorrVals=zeros(nCorsPerAv,1);
rng('default');
rng(42);
tic
for i=1:nCorsPerAv
    SelDat=ResponseSpeedNormVect(randperm(nSjs),:);
    MeanOddVect=mean(SelDat(1:2:end,:));
    MeanEvenVect=mean( SelDat(2:2:end,:));
    shCorrVals(i)=corr(MeanOddVect',MeanEvenVect');
end
toc
obsSHcorr=mean(shCorrVals);
shCorrValsPermed=shCorrVals;
meanPermRhos=zeros(nPerms,1);
for perm=1:nPerms
    perm
    tic
    SelDat=ResponseSpeedNorm;
    for sj=1:nSjs
        for im=1:10
            SelDat(sj,:,im)=ResponseSpeedNorm(sj,randperm(1000),im);
        end
    end
    for i=1:nCorsPerAv
        SelDat=SelDat(randperm(48),:);
        MeanOddVect=mean( SelDat(1:2:end,:));
        MeanEvenVect=mean( SelDat(2:2:end,:));
        shCorrValsPermed(i)=corr(MeanOddVect',MeanEvenVect');
    end
    meanPermRhos(perm)=mean(shCorrValsPermed);
    toc
end

if ComputeNow==0
    load ObservedAndPermedAvSHcorrVals_RT.mat % comment this line to recompute
end
figure(444);hold on;
[n, x]=hist(meanPermRhos,80);
hist(meanPermRhos,80);
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [0.5 0.5 0.5];
line([mean( shCorrVals) obsSHcorr], [0 max(n)],'Color',[.8 .1 .1],'LineWidth',3);
xlim([-.03 .09])
ax=gca;
ax.YTick= 0:100:500;
ax.LineWidth=2;
set(gca, 'TickDir', 'out');
obsSHcorr
ProbabilityOfObsSHcorr=mean( meanPermRhos>obsSHcorr)


