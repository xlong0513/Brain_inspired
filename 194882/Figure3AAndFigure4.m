clc
clear all
close all

% *************************************************************************
% Reproduces Figure 3A and Figure 4.
%
% Runs the spiking neural network with all trials.
%
%   Florian Raudies, 09/07/2014, Boston University.
% *************************************************************************

% Set seed for random number generator.
rng(5);

figurePath = './';
LABEL_SIZE = 16;
TITLE_SIZE = 18;

% Simulate the neural network.
nTrial = 130;
[PerCorrect FiringRate FIndex W12perTrial W23perTrial RasterPlot] = ...
    spikingNetworkContextLearning(nTrial);
[nTrial nStim nHippo] = size(FiringRate);


fprintf('Overall percent correct trials: %2.2f.\n',sum(PerCorrect)/nTrial*100);

% *************************************************************************
% Plot the trajectory of percent correct trials.
% *************************************************************************
TrialWindow = repmat(1/30,[30,1]);
PerCorrect  = imfilter(double(PerCorrect),TrialWindow,'same',0)*100;
TrialIndex  = 30:(nTrial-30);
PerCorrect  = PerCorrect(TrialIndex);

figure('Name','Figure 3A', 'NumberTitle','off');
plot(TrialIndex,PerCorrect,'-k','LineWidth',1.5);
xlabel('Sliding 30 Trial Window','FontSize',LABEL_SIZE);
ylabel('Performance (Percent Correct)','FontSize',LABEL_SIZE);
set(gca,'FontSize',LABEL_SIZE);
axis([30 nTrial-30 0 110]);
print('-deps',sprintf('%sFigurePercentCorrect%d.eps',figurePath,1));


% *************************************************************************
% Plot the raster plots of cell firing.
% *************************************************************************
PlotIndex   = [1 3 5 7 2 4 6 8];
PlotName    = {'A','B','C','D'};
CellIndex   = [2 5 6 8];
dt          = 0.5;

for iIndex = 1:length(CellIndex),
    iCell = CellIndex(iIndex);
    figure('Position',[50 50 1200 400],'PaperPosition',[1 1 9 3],...
        'Name',['Figure 4',PlotName{iIndex}],'numbertitle','off'); 
    for iStim = 1:nStim,
        iSlot           = sub2ind([nStim nHippo],iStim,iCell);
        DataTrialSpike  = RasterPlot.getAllEntryForSlot(iSlot);
        DataSpike       = DataTrialSpike(:,3:end);
        [Y,X]           = find(DataSpike);
        nSample         = size(DataSpike,1);
        subplot(2,4,PlotIndex(iStim)); 
            plot(X*dt,Y,'.k');
            xlabel('Time to action (ms)');
            ylabel('Sample');
            axis([0 1250 0 50]);
            box on;
            title(sprintf('%s,spikes:%d,samples:%d',...
                index2label(iStim),sum(sum(DataSpike)),nSample));
    end
    print('-deps',sprintf('%sFigureFiringHippoCell%d.eps',figurePath,iCell));
end


% *************************************************************************
% Calculate selectivity index for position / item-position.
% *************************************************************************
nBin = 3;
opt.nTrial      = nTrial;
opt.nStim       = nStim;
opt.nCell       = nHippo;
opt.nBin        = nBin;
[SIPos SIItem SIContext] = firingRateToSI(FiringRate,opt);
