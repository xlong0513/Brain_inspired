clc
clear all
close all

% *************************************************************************
% This script reproduces Figure 3F of the manuscript.
%   Florian Raudies, 01/30/2014, Boston University.
%   This script will run for about 45 minutes.
% *************************************************************************

LABEL_SIZE = 16;

NHidden     = 10:10:80;
nNHidden    = length(NHidden);
nLayer      = 3;
nRun        = 50;
nBlock      = 200;
Err         = zeros(nRun,nNHidden);
LetterLabel = {'A','B','C','D'};
NumberLabel = {'1','2','3','4'};

tic
for iRun = 1:nRun,
    % Set seed for random number generator to be able to replicate data.
    rng(1+iRun);
    fprintf('Working on run %d of %d.\n',iRun,nRun);
    for iHidden = 1:nNHidden,
        fprintf('Working on hidden number %d of %d.\n',iHidden,nNHidden);
        nHidden = NHidden(iHidden);
        dcl = DoubleContextLearnerDBNaLP(LetterLabel,NumberLabel,NHidden(iHidden),nLayer);
        dcl.learn(nBlock,{'A1','B1'});
        Err(iRun,iHidden) = dcl.testError;
    end
end
toc

nSample = size(Err,1);
sErr    = 1/sqrt(nSample);
id      = dcl.getIdentifier();

figure('Position',[50 50 600 500],'PaperPosition',[2 2 5 4],'Name','F');
bar(NHidden,mean(Err,1),'FaceColor',[0.7 0.7 0.7]); hold on;
errorbar(NHidden,mean(Err,1),sErr*std(Err,0,1),'k.',...
        'LineWidth',1.5); hold off;
xlabel('Number of hidden neurons','FontSize',LABEL_SIZE);
ylabel('Error probability','FontSize',LABEL_SIZE);
title(id,'FontSize',LABEL_SIZE);
set(gca,'FontSize',LABEL_SIZE);
axis([0 90 0 0.5]); axis square;
