% This program models LGN-V1 pathways
% Feedforward and feedback connections
% Separate excitatory (positive) and inhibitory (negative) connections
% Author: Yanbo Lian
% Date: 26/02/2019
% Citation: Lian Y, Grayden DB, Kameneva T, Meffin H and Burkitt AN (2019) 
%           Toward a Biologically Plausible Model of LGN-V1 Pathways Based
%           on Efficient Coding.
%           Front. Neural Circuits 13:13. doi: 10.3389/fncir.2019.00013

clc; close all; clear

%% Load image
load('IMAGES_SparseCoding.mat')

numImages = size(IMAGES_WHITENED,3);
imageSize = size(IMAGES_WHITENED,1);
imgVar = 0.2; % variance of the input image
BUFF = 4; % the margin between the boundry of the image and selected patch

histFlag = 1; % display the history of cells responses
displayEvery = 200; % display plots after some trials
resizeFactor = 3; % higher resolution when displaying images

%% Define hyper parameters
numPretrain = 1e4;
numEpoches = 3e4; % number of epoches
batchSize = 100; % number of natural images in a minibatch
batchSizePretrain = 100; % number of images of white noise in a minibatch

normalizationMethod = 'L2 norm'; 
l1 = 1;
l2 = 1;

lambda = 0.6; % control sparseness; threshold of the F-I curve
aEta = 0.5; % learning rate of connections A1
tau = 12; % ms
dt = 3; % ms
uEta = dt/tau; % updating rate of membrane potentials U
nU = 30; % number of iterations of calculating membrane potentials U
threshType = 'non-negative soft'; % type of thresholding function that computes firing rates of simple cells from membrane potentials

%% Definitions of symbols
sz = 16; L = sz^2; % size of the image patch; L ON units and L OFF units
OC = 256/L; % overcompleteness
M1 = OC *L; % number of simple cells

% feedforward (up) connections between 2L LGN cells and M1 simple cells
aInitialMean=0.5; % for exponential distribution: var = mean ^ 2;
initial='exponential'; 
A_Up_Pos = NormalizeA( exprnd(aInitialMean,[2*L M1]), normalizationMethod, l1 ); % positive connections
A_Up_Neg = NormalizeA( -exprnd(aInitialMean,[2*L M1]), normalizationMethod, l2 ); % negative connections
A_Up = A_Up_Pos + A_Up_Neg; % overall feedforward connections

% feedback (down) connections between 2L LGN cells and M1 simple cells
A_Down_Pos = NormalizeA( exprnd(aInitialMean,[2*L M1]), normalizationMethod, l2 ); % positive connections
A_Down_Neg = NormalizeA( -exprnd(aInitialMean,[2*L M1]), normalizationMethod, l1 ); % negative connections
A_Down = A_Down_Pos + A_Down_Neg; % overall feedback connections

dA_Bound = 0.1; % maximal change of synaptic efficacy 

X_Data = zeros( L, batchSize ); % input image patches
X = zeros( 2*L, batchSize ); % input with ON and OFF channels

U_L = randn( 2*L, batchSize ); % membrane potential of ON-OFF LGN cells
S_L = rand( 2*L, batchSize ); % firing rate of ON-OFF LGN cells

U1 = randn( M1, batchSize ); % membrane potential of simple cells
S1 = rand( M1, batchSize ); % firing rate of simple cells

s_b = 2; % background firing rate that gives an offset of the reconstruction error
s1Max = 100; % maximum firing rate of simple cells
sL_Max = 100; % maximum firing rate of LGN cells

errorA_UpDown = ones( 1, 1+numEpoches+numPretrain ); % difference between A_Up and A_Down during learning
errorA_UpPosDownPos = ones( 1, 1+numEpoches+numPretrain ); % difference between A_Up_Pos and A_Down_Pos during learning
errorA_UpNegDownNeg = ones( 1, 1+numEpoches+numPretrain ); % difference between A_Up_Neg and A_Down_Neg during learning
errorA_UpDown(1) = sum ( ( A_Up(:) + A_Down(:) ).^2 ); % initial difference
errorA_UpPosDownPos(1) = sum ( ( A_Up_Pos(:) + A_Down_Neg(:) ).^2 );
errorA_UpNegDownNeg(1) = sum ( ( A_Up_Neg(:) + A_Down_Pos(:) ).^2 );

%% Display A and S

% Display the connections from ON and OFF LGN cells to simple cells
figure(1);
subplot(231); DisplayA( 'ON', A_Up_Pos, resizeFactor ); title('A^{+}_{ON,Up}');
subplot(232); DisplayA( 'ON', A_Up_Neg, resizeFactor ); title('A^{-}_{ON,Up}');
subplot(233); DisplayA( 'ON', A_Up, resizeFactor ); title('A_{ON,Up}');
subplot(234); DisplayA( 'OFF', A_Up_Pos, resizeFactor ); title('A^{+}_{OFF,Up}');
subplot(235); DisplayA( 'OFF', A_Up_Neg, resizeFactor ); title('A^{-}_{OFF,Up}');
subplot(236); DisplayA( 'OFF', A_Up, resizeFactor ); title('A_{OFF,Up}');
colormap(Green2Magenta(64));

% Display the overall receptive fields of simple cells: Aon - Aoff
figure(2);
DisplayA( 'ONOFF', A_Up, resizeFactor); title('RFs: A_{ON,Up}-A_{OFF,Up}');
colormap(scm(256));

% Display the firing rates of LGN cells and simple cells
figure(3);
subplot(211); stem(S_L); title(['S_L: LGN cell responses of ' num2str(batchSize) 'patches']);
xlabel('LGN cells'); ylabel('firing rates');
subplot(212); stem(S1); title(['S1: simple cell responses of ' num2str(batchSize) 'patches']);
xlabel('simple cells'); ylabel('firing rates');

%% Pre-train the model using white noise to make sure A_Up converges to A_Down
X_DataPretrain = zeros( L, batchSizePretrain ); % input image patches
X_Pretrain = zeros( 2*L, batchSizePretrain ); % input with ON and OFF channels

U_L_Pretrain = randn( 2*L, batchSizePretrain ); % membrane potential of ON-OFF LGN cells
S_L_Pretrain = rand( 2*L, batchSizePretrain ); % firing rate of ON-OFF LGN cells

U1Pretrain = randn( M1, batchSizePretrain ); % membrane potential of simple cells
S1Pretrain = rand( M1, batchSizePretrain ); % firing rate of simple cells

for iPretrain = 1 : numPretrain
    
    % Generate white noise input with the variance of 'imgVar'
    X_DataPretrain = sqrt(imgVar) * randn(L, batchSizePretrain);
    
    % ON and OFF LGN input
    X_Pretrain( 1:L, : ) = max( X_DataPretrain, 0 );
    X_Pretrain( L+1:2*L, : ) = -min( X_DataPretrain, 0 );
    
    % Compute S and U for LGN and simple cells using previous values
    [ S1Pretrain, U1Pretrain, S_L_Pretrain, U_L_Pretrain ] = ...
        Compute_S_U_LGN_V1_UpDown( S1Pretrain, U1Pretrain, S_L_Pretrain, U_L_Pretrain,...
            X_Pretrain, A_Up, A_Down, lambda, s_b, uEta, nU, threshType, s1Max, sL_Max);

    % Update up and down connections A1
    dA = aEta * ( S_L_Pretrain - s_b ) * S1Pretrain' / batchSizePretrain; % learning rule
    dA = max( min(dA, dA_Bound), -dA_Bound ); % keep the updated amount bounded
    
    A_Up_Pos = max( A_Up_Pos + 1*dA, 0 );
    A_Up_Neg = min( A_Up_Neg + 1*dA, 0 ); % -A_Up_Neg = max( -A_Up_Neg - dA, 0 );
    
    A_Up_Pos = NormalizeA( A_Up_Pos, normalizationMethod, l1 ); % positive connections
    A_Up_Neg = NormalizeA( A_Up_Neg, normalizationMethod, l2 ); % negative connections
    
    A_Down_Pos = max( A_Down_Pos - 1*dA, 0 ); 
    A_Down_Neg = min( A_Down_Neg - 1*dA, 0 ); % -A_Down_Neg = max( -A_Down_Neg - dA, 0 );
    
    A_Down_Pos = NormalizeA( A_Down_Pos, normalizationMethod, l2 ); % positive connections
    A_Down_Neg = NormalizeA( A_Down_Neg, normalizationMethod, l1 ); % negative connections
    
    A_Up = A_Up_Pos + A_Up_Neg; % overall feedforward connections
    A_Down = A_Down_Pos + A_Down_Neg; % overall feedback connections
    
    max( dA(:) )
    min( dA(:) )

    % Display A and S
    if ( mod(iPretrain,displayEvery) == 0 )
        figure(1); % Display the connections from ON and OFF LGN cells to simple cells
        subplot(231); I_A_ON_Up_Pos = DisplayA( 'ON', A_Up_Pos, resizeFactor ); title('A^{+}_{ON,Up}');
        subplot(232); I_A_ON_Up_Neg = DisplayA( 'ON', A_Up_Neg, resizeFactor ); title('A^{-}_{ON,Up}');
        subplot(233); I_A_ON_Up = DisplayA( 'ON', A_Up, resizeFactor ); title('A_{ON,Up}');
        subplot(234); I_A_OFF_Up_Pos = DisplayA( 'OFF', A_Up_Pos, resizeFactor ); title('A^{+}_{OFF,Up}');
        subplot(235); I_A_OFF_Up_Neg = DisplayA( 'OFF', A_Up_Neg, resizeFactor ); title('A^{-}_{OFF,Up}');
        subplot(236); I_A_OFF_Up = DisplayA( 'OFF', A_Up, resizeFactor ); title('A_{OFF,Up}');
        colormap(Green2Magenta(64));
        
        figure(2); % Display the overall receptive fields of simple cells: Aon - Aoff
        DisplayA( 'ONOFF', A_Up, resizeFactor); title('RFs: A_{ON,Up}-A_{OFF,Up}');
        colormap(scm(256));
        
        figure(3); % Display the firing rates of LGN cells and simple cells
        subplot(211); stem(S_L_Pretrain); title(['S_L: LGN cell responses of ' num2str(batchSize) 'patches']);
        xlabel('LGN cells'); ylabel('firing rates');
        subplot(212); stem(S1Pretrain); title(['S1: simple cell responses of ' num2str(batchSize) 'patches']);
        xlabel('simple cells'); ylabel('firing rates');
    end
    
    % Compute the difference between up and down connections
    errorA_UpDown(1+iPretrain) = sum ( ( A_Up(:) + A_Down(:) ).^2 );
    errorA_UpPosDownPos(1+iPretrain) = sum ( ( A_Up_Pos(:) + A_Down_Neg(:) ).^2 );
    errorA_UpNegDownNeg(1+iPretrain) = sum ( ( A_Up_Neg(:) + A_Down_Pos(:) ).^2 );
    
    % print current status of learning
    fprintf('Pretraining %6d: ||Aup-Adown||^2: %4.4f\n',...
                iPretrain, errorA_UpDown(1+iPretrain));
% 	pause
end

        
%% train the model using whitened natural images
for iEpoch = 1 : numEpoches
    
    % adjust the learning rate
    if iEpoch > 1e4
        aEta = 0.2;
    end
    
    if iEpoch > 2e4
        aEta = 0.1;
    end
    
    % Choose an image at random out of 10 images in the dataset
    iImage = ceil( numImages * rand );
    thisImage = IMAGES_WHITENED(:,:,iImage);
    
    % extract image patches at random from this image to make data vector
    for iBatch = 1 : batchSize
        r = BUFF + ceil((imageSize-sz-2*BUFF)*rand); % select y coordinate
        c = BUFF + ceil((imageSize-sz-2*BUFF)*rand); % select x coordinate
        X_Data( : , iBatch ) = reshape( thisImage(r:r+sz-1,c:c+sz-1), L, 1 );
    end
    
    % ON and OFF LGN input
    X_ON = max( X_Data, 0 );
    X_OFF = -min( X_Data, 0 );
    
    X( 1:L, : ) = X_ON;
    X( L+1:2*L, : ) = X_OFF;
    
    % Compute S and U for LGN and simple cells using previous values
    [ S1, U1, S_L, U_L, S1_hist] = Compute_S_U_LGN_V1_UpDown( S1, U1, S_L, U_L,...
        1*X, A_Up, A_Down, lambda, s_b, uEta, nU, threshType, s1Max, sL_Max, histFlag);

    % Update up and down connections A1
    dA = aEta * ( S_L - s_b ) * S1' / batchSize; % learning rule
    dA = max( min(dA, dA_Bound), -dA_Bound ); % keep the updated amount bounded
    
    A_Up_Pos = max( A_Up_Pos + 1*dA, 0 );
    A_Up_Neg = min( A_Up_Neg + 1*dA, 0 ); % -A_Up_Neg = max( -A_Up_Neg - dA, 0 );
    A_Up_Pos = NormalizeA( A_Up_Pos, normalizationMethod, l1 ); % positive connections
    A_Up_Neg = NormalizeA( A_Up_Neg, normalizationMethod, l2 ); % negative connections
    
    A_Down_Pos = max( A_Down_Pos - 1*dA, 0 ); 
    A_Down_Neg = min( A_Down_Neg - 1*dA, 0 );
    A_Down_Pos = NormalizeA( A_Down_Pos, normalizationMethod, l2 ); % positive connections
    A_Down_Neg = NormalizeA( A_Down_Neg, normalizationMethod, l1 ); % negative connections
    
    A_Up = A_Up_Pos + A_Up_Neg; % overall feedforward connections
    A_Down = A_Down_Pos + A_Down_Neg; % overall feedback connections
    
    max( dA(:) )
    min( dA(:) )

    % Display A and S
    if ( mod(iEpoch,displayEvery) == 0 )
        figure(1); % Display the connections from ON and OFF LGN cells to simple cells
        subplot(231); DisplayA( 'ON', A_Up_Pos, resizeFactor ); title('A^{+}_{ON,Up}');
        subplot(232); DisplayA( 'ON', A_Up_Neg, resizeFactor ); title('A^{-}_{ON,Up}');
        subplot(233); DisplayA( 'ON', A_Up, resizeFactor ); title('A_{ON,Up}');
        subplot(234); DisplayA( 'OFF', A_Up_Pos, resizeFactor ); title('A^{+}_{OFF,Up}');
        subplot(235); DisplayA( 'OFF', A_Up_Neg, resizeFactor ); title('A^{-}_{OFF,Up}');
        subplot(236); DisplayA( 'OFF', A_Up, resizeFactor ); title('A_{OFF,Up}');
        colormap(Green2Magenta(64));
        
        figure(2); % Display the overall receptive fields of simple cells: Aon - Aoff
        DisplayA( 'ONOFF', A_Up, resizeFactor); title('RFs: A_{ON,Up}-A_{OFF,Up}');
        colormap(scm(256));
        
        figure(3); % Display the firing rates of LGN cells and simple cells
        subplot(211); stem(S_L); title(['S_L: LGN cell responses of ' num2str(batchSize) 'patches']);
        xlabel('LGN cells'); ylabel('firing rates');
        subplot(212); stem(S1); title(['S1: simple cell responses of ' num2str(batchSize) 'patches']);
        xlabel('simple cells'); ylabel('firing rates');
        
        % Display the trajectory of simple cells responses
        if histFlag == 1
            figure(4);
            plot(S1_hist);title('Trajectory of simple cells')
        end
    end
    
    % Compute the difference between up and down connections
    errorA_UpDown(1+numPretrain+iEpoch) = sum ( ( A_Up(:) + A_Down(:) ).^2 );
    errorA_UpPosDownPos(1+numPretrain+iEpoch) = sum ( ( A_Up_Pos(:) + A_Down_Neg(:) ).^2 );
    errorA_UpNegDownNeg(1+numPretrain+iEpoch) = sum ( ( A_Up_Neg(:) + A_Down_Pos(:) ).^2 );
    
    % print current status of learning
    fprintf('Iteration %6d: ||Aup-Adown||^2: %4.4f\n',...
                iEpoch,errorA_UpDown(1+numPretrain+iEpoch));
end
