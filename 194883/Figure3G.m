clc
clear all
close all

% *************************************************************************
% This script reproduces Figure 3G of the manuscript.
%   Florian Raudies, 01/30/2014, Boston University.
%   This script will run for about 35 minutes.
% *************************************************************************

LABEL_SIZE = 16;

NLayer      = 1:6;
nNLayer     = length(NLayer);
nHidden     = 40;
nRun        = 50;
nBlock      = 200;
Err         = zeros(nRun,nNLayer);
LetterLabel = {'A','B','C','D'};
NumberLabel = {'1','2','3','4'};

tic
for iRun = 1:nRun,
    % Set seed for random number generator to be able to replicate data.
    rng(1+iRun);
    fprintf('Working on run %d of %d.\n',iRun,nRun);
    for iLayer = 1:nNLayer,
        fprintf('Working on block number %d of %d.\n',iLayer,nNLayer);
        dcl = DoubleContextLearnerDBNaLP(LetterLabel,NumberLabel,nHidden,NLayer(iLayer));
        dcl.learn(nBlock,{'A1','B1'});
        Err(iRun,iLayer) = dcl.testError;
    end
end
toc

nSample = size(Err,1);
sErr    = 1/sqrt(nSample);
id      = dcl.getIdentifier();

figure('Position',[50 50 1200 500],'PaperPosition',[2 2 10 4],'Name','3G');
bar(NLayer,mean(Err,1),'FaceColor',[0.7 0.7 0.7]); hold on;
errorbar(NLayer,mean(Err,1),sErr*std(Err,0,1),'k.',...
        'LineWidth',1.5); hold off;
xlabel('Number of layers','FontSize',LABEL_SIZE);
ylabel('Error probability','FontSize',LABEL_SIZE);
title(id,'FontSize',LABEL_SIZE);
set(gca,'FontSize',LABEL_SIZE);
axis([0 7 0 0.5]); axis square;
