clc
clear all
close all

% *************************************************************************
% This script reproduces Figure 3C of the manuscript.
%   Florian Raudies, 01/30/2014, Boston University.
%   This script will run for about 50 minutes.
% *************************************************************************

LABEL_SIZE = 16;

nLayer      = 3;
nHidden     = 40;
nRun        = 50;
NBlock      = [100 150 200 300 400 800];
nNBlock     = length(NBlock);
dcl         = DoubleContextLearnerDBNaLP({'A','B','C','D'},...
                                         {'1','2','3','4'},nHidden,nLayer);

Err  = zeros(nRun,nNBlock);
tic
for iRun = 1:nRun,
    % Set seed for random number generator to be able to replicate data.
    rng(1+iRun);
    fprintf('Working on run %d of %d.\n',iRun,nRun);
    for iBlock = 1:nNBlock,
        fprintf('Working on block number %d of %d.\n',iBlock,nNBlock);
        nBlock = NBlock(iBlock);
        dcl.learn(nBlock,{'A1','B1'});
        Err(iRun,iBlock) = dcl.testError;
    end
end
toc

nSample = size(Err,1);
sErr    = 1/sqrt(nSample);
id      = dcl.getIdentifier();

figure('Position',[50 50 600 500],'PaperPosition',[2 2 5 4],'Name','3C');
bar(1:nNBlock,mean(Err,1),'FaceColor',[0.7 0.7 0.7]); hold on;
errorbar(1:nNBlock,mean(Err,1),sErr*std(Err,0,1),'k.',...
        'LineWidth',1.5); hold off;
xlabel('Number of blocks','FontSize',LABEL_SIZE);
ylabel('Error probability','FontSize',LABEL_SIZE);
title(id,'FontSize',LABEL_SIZE);
set(gca,'XTickLabel',num2cellstr(NBlock),'FontSize',LABEL_SIZE);
axis([0 nNBlock+1 0 0.55]); axis square;



