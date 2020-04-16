clc
clear all
close all

% *************************************************************************
% Reproduces Figure 3B and Figure 5.
%
% Runs the spiking neural network with all trials multiple times to collect
% the statistics for selectivity indices. Each run has a different 
% initialization of the random number generator, which results in different 
% weight initlization and noise flucations for the membrane potential.
%
% Use the switch  'SIMULATE' to re-run the simulation, which will take
% about 6-7 hours. For convinience we provide the simulation data.
%
%   Florian Raudies, 09/07/2014, Boston University.
% *************************************************************************

LABEL_SIZE  = 16;
TITLE_SIZE  = 18;
figurePath  = './';
SIMULATE    = 0;

nTrial              = 130;
nRun                = 100;
nBlock              = 30;
PerCorrectPerTrial  = zeros(nRun,nTrial);
FiringRatePerTrial  = cell(nRun,1);
FIndexPerTrial      = cell(nRun,1);
nIn                 = 6;
nHippo              = 8;
nOut                = 2;
nTimeTrial          = 801;
WeightLayer1To2     = zeros(nRun,nTrial,nIn,nHippo);
WeightLayer2To3     = zeros(nRun,nTrial,nHippo,nOut);
fileName            = sprintf('NetworkSimulation%dRuns',nRun);


if SIMULATE,
    
    for iRun = 1:nRun,
        fprintf('Run %d.\n',iRun);
        rng(1+iRun);
        [PerCorrect FiringRate FIndex W12perTrial W23perTrial] ...
            = spikingNetworkContextLearning(nTrial);
        PerCorrectPerTrial(iRun,:)  = PerCorrect;
        FiringRatePerTrial{iRun}    = FiringRate;
        FIndexPerTrial{iRun}        = FIndex;
        WeightLayer1To2(iRun,:,:,:) = W12perTrial;
        WeightLayer2To3(iRun,:,:,:) = W23perTrial;
        fprintf('Overall percent correct trials: %2.2f.\n',...
            sum(PerCorrect)/nTrial*100);
    end

    save(fileName, 'PerCorrectPerTrial', 'FiringRatePerTrial', ...
        'FIndexPerTrial', 'WeightLayer1To2', 'WeightLayer2To3');
    
else
    
    load(fileName);
    
end

% *************************************************************************
% Figure for correct detection graph.
% *************************************************************************
TrialWindow         = repmat(1/30, [1 30]);
PerCorrectPerTrial  = imfilter(double(PerCorrectPerTrial),...
                               TrialWindow,'same',0)*100;
TrialIndex          = 30:100;
PerCorrectPerTrial  = PerCorrectPerTrial(:,TrialIndex);


figure('Name','Figure 3B', 'NumberTitle','off', 'Position',[50 50 800 500]);
errorarea(TrialIndex,mean(PerCorrectPerTrial),...
                     std(PerCorrectPerTrial),[.85 .85 .85],'k');
hold on
plot([30 100],[100 100],'--k');
hold off
axis([30 100 0 110]);
xlabel('Sliding 30 Trial Window','FontSize',LABEL_SIZE);
ylabel('Performance (Percent Correct)','FontSize',LABEL_SIZE);
set(gca,'FontSize',LABEL_SIZE);
title(sprintf('N=%d',nRun),'FontSize',TITLE_SIZE);
print('-deps',sprintf('%sFigurePercentCorrect%dRuns.eps',figurePath,nRun));

% *************************************************************************
% Figures for selectivity indices.
% *************************************************************************
nBin = 4;
[nTrial nStim nHippo]   = size(FiringRatePerTrial{1});
SIPosPerTrial           = zeros(nRun,nBin,nHippo);
SIItemPerTrial          = zeros(nRun,nBin,nHippo);
SIContextPerTrial       = zeros(nRun,nBin,nHippo);
opt.nBin                = nBin;
opt.nTrial              = nTrial;
opt.nCell               = nHippo;
opt.nStim               = nStim;
for iRun = 1:nRun,
    FiringRate                          = FiringRatePerTrial{iRun};
    [SIPos SIItem SIContext]            = firingRateToSI(FiringRate,opt);
    FIndex                              = FIndexPerTrial{iRun};
    SIPosPerTrial(iRun,:,:)             = SIPos;
    SIItemPerTrial(iRun,:,:)            = SIItem;
    SIContextPerTrial(iRun,:,:)         = SIContext;
    SIPosPerTrial(iRun,:,~FIndex)       = NaN;
    SIItemPerTrial(iRun,:,~FIndex)      = NaN;
    SIContextPerTrial(iRun,:,~FIndex)   = NaN;
end
SIPosPerTrial(SIPosPerTrial==0)         = NaN;
SIItemPerTrial(SIItemPerTrial==0)       = NaN;
SIContextPerTrial(SIContextPerTrial==0) = NaN;
MeanSIPos       = squeeze(meanWoutNaN(SIPosPerTrial,1));
SemSIPos        = squeeze(semWoutNaN(SIPosPerTrial,1));
MeanSIItem      = squeeze(meanWoutNaN(SIItemPerTrial,1));
SemSIItem       = squeeze(semWoutNaN(SIItemPerTrial,1));
MeanSIContext   = squeeze(meanWoutNaN(SIContextPerTrial,1));
SemSIContext    = squeeze(semWoutNaN(SIContextPerTrial,1));
XNAME = {'1st','2nd','3rd','4th'};

% *************************************************************************
% Selectivity index for place.
% *************************************************************************
figure('Name','Figur 5A', 'NumberTitle','off', ...
       'Position',[50 50 1200 600],'PaperPosition',[1 1 14 8]);
for iHippo = 1:nHippo,
    subplot(2,nHippo/2,iHippo);
        bar(MeanSIPos(:,iHippo),'EdgeColor',[0 0 0],'FaceColor',[0.7 0.7 0.7]);
        hold on;
        errorbar(1:4,MeanSIPos(:,iHippo),SemSIPos(:,iHippo),'k.','LineWidth',1.5);
        plot([1 4],[1 1],'--k','LineWidth',1.5);
        hold off;
        ylim([0 1.1]);
        set(gca,'XTick',[1 2 3 4],'XTickLabel',XNAME);
        ylabel('SI place','FontSize',LABEL_SIZE);
        set(gca,'FontSize',LABEL_SIZE);
        title(sprintf('Cell %d',iHippo),'FontSize',TITLE_SIZE);
end
print('-depsc',sprintf('%sFigureSIPlace.eps',figurePath));


% *************************************************************************
% Paired t-test. The significance interval is 0.05 or 5%.
% *************************************************************************
fprintf('ttest2 for first and last 30 trials based on SI for place.\n');
for iHippo = 1:nHippo,
    hypo = ttest2(SIPosPerTrial(:,1,iHippo),SIPosPerTrial(:,4,iHippo),.05);
    fprintf('Hippocampal cell %d is significantly different %d.\n',iHippo,hypo);
end


% *************************************************************************
% Selectivity index for item / place.
% *************************************************************************
figure('Name','Figure 5B', 'NumberTitle','off', ...
       'Position',[50 50 1200 600],'PaperPosition',[1 1 14 8]); 
for iHippo = 1:nHippo,
    subplot(2,nHippo/2,iHippo);
        bar(MeanSIItem(:,iHippo),'EdgeColor',[0 0 0],'FaceColor',[0.7 0.7 0.7]);
        hold on;
        errorbar(1:4,MeanSIItem(:,iHippo),SemSIItem(:,iHippo),'k.','LineWidth',1.5);
        plot([1 4],[1 1],'--k','LineWidth',1.5);
        hold off;
        ylim([0 1.1]);
        set(gca,'XTick',[1 2 3 4],'XTickLabel',XNAME);
        ylabel('SI item','FontSize',LABEL_SIZE);
        set(gca,'FontSize',LABEL_SIZE);
        title(sprintf('Cell %d',iHippo),'FontSize',TITLE_SIZE);
end
print('-depsc',sprintf('%sFigureSIItem.eps',figurePath));


% *************************************************************************
% Paired t-test. The significance interval is 0.01 or 1%.
% *************************************************************************
fprintf('ttest2 for first and last 30 trials based on SI for item.\n');
for iHippo = 1:nHippo,
    hypo = ttest2(SIItemPerTrial(:,1,iHippo),SIItemPerTrial(:,4,iHippo),.01);
    fprintf('Hippocampal cell %d is significantly different %d.\n',iHippo,hypo);
end


% *************************************************************************
% Selectivity index for context.
% *************************************************************************
figure('Name','Figure 5C','NumberTitle','off', ...
       'Position',[50 50 1200 600],'PaperPosition',[1 1 14 8]);
for iHippo = 1:nHippo,
    subplot(2,nHippo/2,iHippo);
        bar(MeanSIContext(:,iHippo),'EdgeColor',[0 0 0],'FaceColor',[0.7 0.7 0.7]);
        hold on;
        errorbar(1:4,MeanSIContext(:,iHippo),SemSIContext(:,iHippo),'k.','LineWidth',1.5);
        plot([1 4],[1 1],'--k','LineWidth',1.5);
        hold off;
        ylim([0 1.1]);
        set(gca,'XTick',[1 2 3 4],'XTickLabel',XNAME);
        ylabel('SI context','FontSize',LABEL_SIZE);
        set(gca,'FontSize',LABEL_SIZE);
        title(sprintf('Cell %d',iHippo),'FontSize',TITLE_SIZE);
end
print('-depsc',sprintf('%sFigureSIContext.eps',figurePath));

% *************************************************************************
% Paired t-test. The significance interval is 0.01 or 1%.
% *************************************************************************
fprintf('ttest2 for first and last 30 trials based on SI for context.\n');
for iHippo = 1:nHippo,
    hypo = ttest2(SIContextPerTrial(:,1,iHippo),SIContextPerTrial(:,4,iHippo),.01);
    fprintf('Hippocampal cell %d is significantly different %d.\n',iHippo,hypo);
end


% *************************************************************************
% Compute the binariness for each of the blocks using the weights that 
% connect layer 1 to layer 2.
%   Florian Raudies, 07/23/2014, Boston University.
% *************************************************************************
Binariness12    = binariness(WeightLayer1To2);
Binariness23    = binariness(WeightLayer2To3);
FIndexPerTrial  = cell2mat(FIndexPerTrial);
FIndexPerTrial  = permute(repmat(FIndexPerTrial,[1 1 nTrial nIn]),[1 3 4 2]);
Binariness12(~FIndexPerTrial) = NaN;
Binariness12 = squeeze(meanWoutNaN(Binariness12,3));
Binariness12PerBin = zeros(nRun,nBin,nHippo);
for iBin = 1:nBin,
    Binariness12PerBin(:,iBin,:) = meanWoutNaN(...
        Binariness12(:,(iBin-1)*nBlock+(1:nBlock),:),2);
end
MeanBinariness = squeeze(meanWoutNaN(Binariness12PerBin,1));
SemBinariness = squeeze(semWoutNan(Binariness12PerBin,1));


figure('Name','Figure 5D','NumberTitle','off',...
       'Position',[50 50 1200 600],'PaperPosition',[1 1 14 8]);
for iHippo = 1:nHippo,
    subplot(2,nHippo/2,iHippo);
        bar(MeanBinariness(:,iHippo),'EdgeColor',[0 0 0],'FaceColor',[0.7 0.7 0.7]);
        hold on;
        errorbar(1:4,MeanBinariness(:,iHippo),SemBinariness(:,iHippo),'k.','LineWidth',1.5);
        plot([1 4],[1 1],'--k','LineWidth',1.5);
        hold off;
        ylim([0 1.1]);
        set(gca,'XTick',[1 2 3 4],'XTickLabel',XNAME);
        ylabel('Binariness','FontSize',LABEL_SIZE);
        set(gca,'FontSize',LABEL_SIZE);
        title(sprintf('Cell %d',iHippo),'FontSize',TITLE_SIZE);
end
print('-depsc',sprintf('%sFigureSIBinariness.eps',figurePath));


% *************************************************************************
% Paired t-test for binariness. The significance interval is 0.01 or 1%.
% *************************************************************************
fprintf('ttest2 for first and last 30 trials based on binariness.\n');
for iHippo = 1:nHippo,
    hypo = ttest2(Binariness12PerBin(:,1,iHippo),Binariness12PerBin(:,4,iHippo),.01);
    fprintf('Hippocampal cell %d is significantly different %d.\n',iHippo,hypo);
end
