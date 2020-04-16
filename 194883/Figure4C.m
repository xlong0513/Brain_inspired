clc
clear all
close all

% *************************************************************************
% This script reproduces Figure 4C of the manuscript.
%   Florian Raudies, 01/30/2014, Boston University.
%   This script will run for about 1,100 minutes or 16 hours.
% *************************************************************************

LABEL_SIZE = 16;

ExcludeList = {{'A1','B1'},{'A1','A2'},{'A1','A4'},...
               {'A1','B4'},{'A1','C1'},{'A1','C4'},...
               {'A1','A2','A3'},{'A1','B1','C1'},{'A1','A2','B2'},...
               {'A1','C2','D2'},{'A1','B2','D2'},{'A1','A2','A4'},...
               {'A1','B1','C1','D1'},{'A1','A2','A3','A4'},...
               {'A1','B1','A2','B2'},{'A1','B1','A4','B4'},...
               {'A1','C1','A2','C2'},{'A1','D1','A4','D4'}};
nExclude    = size(ExcludeList,2);
nHidden     = 40;
nRun        = 50;
nBlock      = 200;
LetterLabel = {'A','B','C','D'};
NumberLabel = {'1','2','3','4'};
dcl         = DoubleContextLearnerMLP(LetterLabel,NumberLabel,nHidden);
Err         = zeros(nRun,nExclude);

tic
for iRun = 1:nRun,
    % Set seed for random number generator to be able to replicate data.
    rng(1+iRun);
    fprintf('Working on run %d of %d.\n',iRun,nRun);
    for iExclude = 1:nExclude,
        fprintf('Working on exlcude case %d of %d.\n',iExclude,nExclude);
        dcl.learn(nBlock,ExcludeList{iExclude});
        Err(iRun,iExclude) = dcl.testError;
    end
end
toc


LabelList = cell(nExclude,1);
for iExclude = 1:nExclude,
    ExcludeItem = ExcludeList{iExclude};
    nItem = length(ExcludeItem);
    LabelList{iExclude} = '';
    for iItem = 1:nItem,
        LabelList{iExclude} = [LabelList{iExclude}, ExcludeItem{iItem}];
    end
end


nSample = size(Err,1);
sErr    = 1/sqrt(nSample);
id      = dcl.getIdentifier();

figure('Position',[50 50 600 600],'PaperPosition',[2 2 5 6],'Name','4C');
bar(1:nExclude,mean(Err,1),'FaceColor',[0.7 0.7 0.7]); hold on;
errorbar(1:nExclude,mean(Err,1),sErr*std(Err,0,1),'k.',...
        'LineWidth',1.5); hold off;
xlabel('Excluded from training','FontSize',LABEL_SIZE);
ylabel('Error probability','FontSize',LABEL_SIZE);
title(id,'FontSize',LABEL_SIZE);
set(gca,'XTick',1:nExclude,'XTickLabel',LabelList,'FontSize',LABEL_SIZE);
rotateXLabels(gca, 90);
axis([0 nExclude+1 0 0.6]); axis square;

