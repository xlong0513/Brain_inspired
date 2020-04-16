classdef RestrictedBoltzmannMachine < handle
    % Restricted Boltzmann Machine (RBM) with visible and hidden layer and
    % weights in between.
    % Uses contrastive divergence to adapt weights and thresholds. If this
    % RBM is used as last layer with labels use the method
    % fitContrastiveDivergence to adapt weights and thresholds.
    %
    %   Florian Raudies, 01/08/2014, Boston University.
    %
    %   This is a re-implementation of Andrej Karpathy's code which was 
    %   based on Kevin Swersky and Ruslan Salakhutdinov's code.
    properties (SetAccess = private)
        W   % Weights between visible and hidden layer.
        C   % Biases of the visible layer.
        B   % Biases of the hidden layer.
        E   % Training error for each epoch.
        T   % Top layer activity.
        Wc  % Weights for classes.
        Cc  % Biases for classes.
    end
    properties
        nMaxEpoch   % Maximum number of epochs.
        nAvg        % Number of averaging epochs.
        nBatchSize  % Number of trials per batch.
        eta         % Learning rate.
        momentum    % Momentum term.
        penalty     % Weight decay factor.
        verbose     % Print info.
        anneal      % Use annealing.
        blockTrain  % Train with ordered block.
    end
    methods
        % Constructor
        function obj = RestrictedBoltzmannMachine()
            obj.nMaxEpoch   = 50;
            obj.nAvg        = 5;
            obj.nBatchSize  = 100;
            obj.eta         = 0.2;
            obj.momentum    = 0.5;
            obj.penalty     = 2e-4;
            obj.verbose     = 0;
            obj.anneal      = 0;
            obj.blockTrain  = 0;
        end
        % Train using the contrastive convergence (CC) algorithm.
        function obj = trainContrastiveConvergence(obj, Data, nHidden)
            % Split data into batches.
            nData       = size(Data,1);
            nDim        = size(Data,2);
            nBatch      = ceil(nData/obj.nBatchSize);
            BatchData   = cell(nBatch,1);
            % Bring data into random order.
            % Note that this shuffling schema for blocks is important. A
            % simple shuffle applying randperm to the entire data set and
            % piecing it into nBatchSize junks does not work. FR 01/07/2014.
            if obj.blockTrain
                for iBatch = 1:nBatch-1
                    BatchData{iBatch} = Data(...
                        (iBatch-1)*obj.nBatchSize+(1:obj.nBatchSize),:);
                end
                BatchData{nBatch} = Data( ...
                    (nBatch-1)*obj.nBatchSize+1:nData,:);
            else
                Shuffle     = repmat(1:nBatch, [1 obj.nBatchSize]);
                Shuffle     = Shuffle(randperm(nData));
                for iBatch = 1:nBatch,
                    BatchData{iBatch} = Data(Shuffle==iBatch,:);
                end
            end
            % Initalize the variables.
            obj.W    = 0.01 * randn(nDim, nHidden);
            obj.C    = repmat( log(.25/.75), [1 nDim]);
            obj.B    = zeros(1, nHidden);
            Winc     = zeros(nDim,nHidden);
            Binc     = zeros(1,nHidden);
            Cinc     = zeros(1,nDim);
            Wavg     = obj.W;
            Bavg     = obj.B;
            Cavg     = obj.C;
            p        = obj.penalty;
            iAvg     = 1;
            obj.E    = zeros(obj.nMaxEpoch,1);
            nAvgStart = obj.nMaxEpoch - obj.nAvg;
            for iEpoch = 1:obj.nMaxEpoch,
                errSum = 0;
                if (obj.anneal)
                    p = obj.penalty - 0.9*iEpoch/obj.nMaxEpoch*obj.penalty;
                end               
                for iBatch = 1:nBatch,
                    DataBatch = BatchData{iBatch};
                    nThisBatchSize = size(DataBatch,1);
                    % Perform bottom-up pass.
                    PosHid = logistic(DataBatch*obj.W ...
                                   + repmat(obj.B,[nThisBatchSize 1]));
                    PosHidStates = PosHid > rand(nThisBatchSize, nHidden);
                    % Negative phase.
                    NegData = logistic(PosHidStates*obj.W' ...
                                     + repmat(obj.C,[nThisBatchSize 1]));
                    NegDataStates = NegData > rand(nThisBatchSize, nDim);
                    % Perform 2nd bottom-up pass.
                    NegHid = logistic(NegDataStates*obj.W ...
                                     + repmat(obj.B,[nThisBatchSize 1]));
%                    HegHidStates = NegHid > rand(obj.nBatchSize, nHidden);
                    % Update weights.
                    DW = DataBatch'*PosHid - NegDataStates'*NegHid;
                    DC = sum(DataBatch) - sum(NegDataStates);
                    DB = sum(PosHid) - sum(NegHid);
                    Winc = obj.momentum * Winc ...
                         + obj.eta*(DW/nThisBatchSize - p*obj.W);
                    Binc = obj.momentum * Binc ...
                         + obj.eta*(DB/nThisBatchSize);
                    Cinc = obj.momentum * Cinc ...
                         + obj.eta*(DC/nThisBatchSize);
                    obj.W = obj.W + Winc;
                    obj.B = obj.B + Binc;
                    obj.C = obj.C + Cinc;
                    % Average weights from the last nMaxEpoch-nAvgStart+1
                    % epochs.
                    if iEpoch > nAvgStart,
                        Wavg = Wavg - (1/iAvg)*(Wavg - obj.W);
                        Cavg = Cavg - (1/iAvg)*(Cavg - obj.C);
                        Bavg = Bavg - (1/iAvg)*(Bavg - obj.B);
                        iAvg = iAvg + 1;
                    else
                        Wavg = obj.W;
                        Bavg = obj.B;
                        Cavg = obj.C;
                    end
                    errSum = errSum + sum(sum( (DataBatch - NegData).^2 ));
                end
                obj.E(iEpoch) = errSum;
                if obj.verbose,
                    fprintf('Reconstruction error in epoch %d is %f.\n',...
                        iEpoch, errSum);
                end
            end
            obj.T = logistic(Data*Wavg + repmat(Bavg,[nData 1]));
            obj.W = Wavg;
            obj.B = Bavg;
            obj.C = Cavg;
        end
        function obj = fitContrastiveConvergence(obj,Data,Label,nHidden)
            % Get label and prepare binary data for labels.
            nData    = size(Data,1);
            nDim     = size(Data,2);
            nClasses = size(Label,2);
            % Split data into batches.
            nBatch      = ceil(nData/obj.nBatchSize);
            BatchData   = cell(nBatch,1);
            BatchLabel  = cell(nBatch,1);
            if obj.blockTrain
                for iBatch = 1:nBatch-1
                    Index = (iBatch-1)*obj.nBatchSize+(1:obj.nBatchSize);
                    BatchData{iBatch} = Data(Index,:);
                    BatchLabel{iBatch} = Label(Index,:);
                end
                Index = (nBatch-1)*obj.nBatchSize+1:nData;
                BatchData{nBatch} = Data(Index,:);
                BatchLabel{nBatch} = Label(Index,:);
            else
                Shuffle     = repmat(1:nBatch, [1 obj.nBatchSize]);
                Shuffle     = Shuffle(randperm(nData));
                for iBatch = 1:nBatch,
                    Index = Shuffle==iBatch;
                    BatchData{iBatch} = Data(Index,:);
                    BatchLabel{iBatch} = Label(Index,:);
                end
            end
            % Initalize the variables.
            obj.W    = 0.01 * randn(nDim, nHidden);
            obj.C    = repmat( log(.25/.75), [1 nDim]);
            obj.B    = zeros(1, nHidden);
            obj.Wc   = 0.01 * randn(nClasses, nHidden);
            obj.Cc   = zeros(1, nClasses);
            Winc     = zeros(nDim,nHidden);
            Binc     = zeros(1,nHidden);
            Cinc     = zeros(1,nDim);
            Wcinc    = zeros(nClasses, nHidden);
            Ccinc    = zeros(1, nClasses);
            Wavg     = obj.W;
            Bavg     = obj.B;
            Cavg     = obj.C;
            Wcavg    = obj.Wc;
            Ccavg    = obj.Cc;
            p        = obj.penalty;
            iAvg     = 1;
            obj.E    = zeros(obj.nMaxEpoch,1);
            nAvgStart = obj.nMaxEpoch - obj.nAvg;
            for iEpoch = 1:obj.nMaxEpoch,
                errSum = 0;
                if (obj.anneal)
                    p = obj.penalty - 0.9*iEpoch/obj.nMaxEpoch*obj.penalty;
                end               
                for iBatch = 1:nBatch,
                    DataBatch = BatchData{iBatch};
                    LabelBatch = BatchLabel{iBatch};
                    nThisBatchSize = size(DataBatch,1);
                    % Perform bottom-up pass.
                    PosHid = logistic(DataBatch*obj.W + LabelBatch*obj.Wc ...
                                   + repmat(obj.B,[nThisBatchSize 1]));
                    PosHidStates = PosHid > rand(nThisBatchSize, nHidden);
                    % Negative phase
                    NegData = logistic(PosHidStates*obj.W' ...
                                     + repmat(obj.C,[nThisBatchSize 1]));
                    NegDataStates = NegData > rand(nThisBatchSize, nDim);
                    NegLabel = RestrictedBoltzmannMachine.softmaxPmtk( ...
                                 PosHidStates*obj.Wc' ...
                               + repmat(obj.Cc,[nThisBatchSize 1]));
                    NegLabelStates = RestrictedBoltzmannMachine.softmaxSample(...
                                 NegLabel);
                    % Perform 2nd bottom-up pass
                    NegHid = logistic(NegDataStates*obj.W + NegLabelStates*obj.Wc + ...
                                     + repmat(obj.B,[nThisBatchSize 1]));
%                    HegHidStates = NegHid > rand(obj.nBatchSize, nHidden);
                    % Update weights
                    DW       = DataBatch'*PosHid - NegDataStates'*NegHid;
                    DC       = sum(DataBatch) - sum(NegDataStates);
                    DB       = sum(PosHid) - sum(NegHid);
                    DWc      = LabelBatch'*PosHid - NegLabelStates'*NegHid;
                    DCc      = sum(LabelBatch) - sum(NegLabelStates);
                    Winc     = obj.momentum * Winc ...
                             + obj.eta*(DW/nThisBatchSize - p*obj.W);
                    Binc     = obj.momentum * Binc ...
                             + obj.eta*(DB/nThisBatchSize);
                    Cinc     = obj.momentum * Cinc ...
                             + obj.eta*(DC/nThisBatchSize);
                    Wcinc    = obj.momentum * Wcinc ...
                             + obj.eta*(DWc/nThisBatchSize - p*obj.Wc);
                    Ccinc    = obj.momentum * Ccinc ...
                             + obj.eta*DCc/nThisBatchSize;
                    obj.W    = obj.W + Winc;
                    obj.B    = obj.B + Binc;
                    obj.C    = obj.C + Cinc;
                    obj.Wc   = obj.Wc + Wcinc;
                    obj.Cc   = obj.Cc + Ccinc;
                    % Average weights from the last nMaxEpoch-nAvgStart+1
                    % epochs.
                    if iEpoch > nAvgStart,
                        Wavg     = Wavg - (1/iAvg)*(Wavg - obj.W);
                        Cavg     = Cavg - (1/iAvg)*(Cavg - obj.C);
                        Bavg     = Bavg - (1/iAvg)*(Bavg - obj.B);
                        Wcavg    = Wcavg - (1/iAvg)*(Wcavg - obj.Wc);
                        Ccavg    = Ccavg - (1/iAvg)*(Ccavg - obj.Cc);
                        iAvg     = iAvg + 1;
                    else
                        Wavg = obj.W;
                        Bavg = obj.B;
                        Cavg = obj.C;
                    end
                    errSum = errSum + sum(sum( (DataBatch - NegData).^2 ));
                end
                obj.E(iEpoch) = errSum;
                if obj.verbose,
                    fprintf('Reconstruction error in epoch %d is %f.\n',...
                        iEpoch, errSum);
                end
            end
            obj.W    = Wavg;
            obj.B    = Bavg;
            obj.C    = Cavg;
            obj.Wc   = Wcavg;
            obj.Cc   = Ccavg;
        end
        % For the last layer.
        function P = predict(obj,Data)
            if isempty(obj.Wc) || isempty(obj.Cc),
                error('Matlab:RestrictedBoltzmannMachine',...
                    'No prediction possible. This is not the output layer\n');
            end
            nClasses = size(obj.Wc,1);
            nData = size(Data,1);
            F = zeros(nData,nClasses);
            for iClasses = 1:nClasses,
                X = zeros(nData,nClasses);
                X(:,iClasses) = 1;
                F(:,iClasses) = repmat(obj.Cc(iClasses),[nData 1]) ...
                    + sum(log(1+exp(Data*obj.W + X*obj.Wc + repmat(obj.B,[nData 1]))),2);
            end
            [~, Index] = max(F,[],2);
            P = zeros(nData,nClasses);
            Index = sub2ind([nData nClasses], (1:nData)',Index);            
            P(Index) = 1;
        end
        function E = trainingError(obj)
            E = obj.E;
        end
        function H = visibleToHidden(obj,V)
            H = logistic(V*obj.W + repmat(obj.B,size(V,1),1));
        end
        function V = hiddenToVisible(obj,H)
            V = logistic(H*obj.W' + repmat(obj.C,size(H,1),1));
        end
        function W = getWeight(obj)
            W = obj.W;
        end
        function Wc = getWeightForClass(obj)
            Wc = obj.Wc;
        end
        function obj = setBlockTrain(obj,flag)
            obj.blockTrain = flag;
        end
    end
    methods (Static = true)
        % Calculate softmax value.
        function X = softmaxPmtk(X)
            X = exp(X);
            X = bsxfun(@rdivide, X, sum(X, 2));
        end
        % One sample per row.
        function S = softmaxSample(P)
            nRow = size(P,1);
            P    = bsxfun(@rdivide, P, sum(P, 2)); % normalize
            Th   = bsxfun(@gt,cumsum(P,2), rand([nRow,1])); % threshold
            S    = double(diff([zeros(nRow,1),Th],1,2)>0); % get first occurence
        end
    end
end