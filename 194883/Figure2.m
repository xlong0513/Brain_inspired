clc
clear all
close all

% *************************************************************************
% This script reproduces Figure 2B-G of the manuscript.
%   Florian Raudies, 01/30/2014, Boston University.
%   This script will run for about 2 minutes.
% *************************************************************************

LABEL_SIZE = 16;

% Set seed for random number generator to be able to replicate data.
rng(2);

LetterLabel = {'A','B','C','D'};
NumberLabel = {'1','2','3','4'};

% Initialize the four quadrant learner.
dcl = DoubleContextLearnerDBNaLP(LetterLabel,NumberLabel,40,3);
dcl.learn(800,{'A1','B1'});
err = dcl.testError;

[A Wc] = dcl.getDBNActivationSortedByWeights();
[nLayer nHidden nState] = size(A);
PlotIndex = [1  5  2  6; 9  13 10 14; 3  7  4  8; 11 15 12 16];
Asum = cumsum(repmat(Wc,[nState 1])'.*squeeze(A(end,:,:)));
Data = sum(repmat(Wc,[nState 1])'.*squeeze(A(end,:,:)));

A = (A-min(A(:)))/(max(A(:))-min(A(:)));
Asum = (Asum-min(Asum(:)))/(max(Asum(:))-min(Asum(:)));

% *************************************************************************
% Figure of activations within the network.
% *************************************************************************
figure('Position',[50 50 1200 600],'PaperPosition',[2 2 12 5],'Name','B');
for iLayer = 1:nLayer,
    for iHidden = 1:16,
        iPlot = sub2ind([6 8],2*(iLayer-1)+ceil(iHidden/8),mod(iHidden-1,8)+1);
        subplot(6,8,iPlot);
            imshow(reshape(A(iLayer,iHidden,PlotIndex),[4 4]),[0 1]);
    end
end

% *************************************************************************
% Figure of activations of the network output.
% *************************************************************************
figure('Position',[50 50 1200 200],'PaperPosition',[2 2 12 2],'Name','C');
for iHidden = 1:16,
    iPlot = sub2ind([2 8],ceil(iHidden/8),mod(iHidden-1,8)+1);
    subplot(2,8,iPlot);
        imshow(reshape(Asum(iHidden,PlotIndex),[4 4]),[0 1]);
end

% *************************************************************************
% Figure of the resclaed 16th sum activation.
% *************************************************************************
figure('Name','D');
imshow(reshape(Asum(16,PlotIndex),[4 4]),[],'InitialMagnification',3*10^3);

Label = dcl.getLabelBlock();
[~, Label] = max(Label,[],2);
Label = Label - 1;
% Determine class labels for display purposes.
CY = Label > 0;
CX = ~CY;
% Retrieve parameters and visualize the result.
w       = dcl.lp.getWeight;
theta   = dcl.lp.getThreshold;
X       = -10:2;
Y       = -w*X+theta;

% *************************************************************************
% Figure of the linear perceptron hyperplane and data.
% *************************************************************************
figure('Name','E');
h1 = plot(Data(CX),0,'k^','MarkerSize',12); hold on;
h2 = plot(Data(CY),0,'ko','MarkerSize',12);
h3 = plot(X,Y,'-k','LineWidth',1.5); hold off;
legend([h1(1) h2(2) h3],'X','Y','Hyperplane','Location','SouthEast');
axis equal; axis([-10 2 -.5 .5]); 
xlabel('activation','FontSize',LABEL_SIZE);
ylabel('auxiliary','FontSize',LABEL_SIZE);
set(gca,'FontSize',LABEL_SIZE);


% *************************************************************************
% Train with all stimuli.
% *************************************************************************
rng(2);
% Initialize the four quadrant learner.
dcl = DoubleContextLearnerDBNaLP(LetterLabel,NumberLabel,40,3);
dcl.learn(800,{});
err = dcl.testError;

[A Wc] = dcl.getDBNActivationSortedByWeights();
[nLayer nHidden nState] = size(A);
Asum = cumsum(repmat(Wc,[nState 1])'.*squeeze(A(end,:,:)));
Data = sum(repmat(Wc,[nState 1])'.*squeeze(A(end,:,:)));

% *************************************************************************
% Figure of the resclaed 16th sum activation.
% *************************************************************************
figure('Name','F');
imshow(reshape(Asum(16,PlotIndex),[4 4]),[],'InitialMagnification',3*10^3);

Label = dcl.getLabelBlock();
[~, Label] = max(Label,[],2);
Label = Label - 1;
% Determine class labels for display purposes.
CY = Label > 0;
CX = ~CY;
% Retrieve parameters and visualize the result.
w       = dcl.lp.getWeight;
theta   = dcl.lp.getThreshold;
X       = -10:2;
Y       = -w*X+theta;

% *************************************************************************
% Figure of the linear perceptron hyperplane and data.
% *************************************************************************
figure('Name','G');
h1 = plot(Data(CX),0,'k^','MarkerSize',12); hold on;
h2 = plot(Data(CY),0,'ko','MarkerSize',12);
h3 = plot(X,Y,'-k','LineWidth',1.5); hold off;
legend([h1(1) h2(2) h3],'X','Y','Hyperplane','Location','SouthEast');
axis equal; axis([-10 2 -.5 .5]); 
xlabel('activation','FontSize',LABEL_SIZE);
ylabel('auxiliary','FontSize',LABEL_SIZE);
set(gca,'FontSize',LABEL_SIZE);

