classdef DoubleContextLearnerDBNaLP < DoubleContextLearner 
    % DoubleContextLearner
    % Uses a the combination of a Deep Belief Network (DBN) and Linear 
    % Perceptron (LP) to learn the double-context task.
    %   Florian Raudies, 01/30/2014, Boston University.
    properties (SetAccess = private)
        dbn         % Deep belief network.
        lp          % Linear perceptron.
    end
    methods
        function obj = DoubleContextLearnerDBNaLP(...
                LetterLabel,NumberLabel,nHidden,nLayer)
            obj = obj@DoubleContextLearner(LetterLabel,NumberLabel);
            obj.dbn = DeepBeliefNetwork(repmat(nHidden,[1 nLayer]));
            obj.lp  = LinearPerceptron;
        end
        function obj = learn(obj,nBlock,ExcludeState)
            [Data Label] = obj.generateData(nBlock, ExcludeState);
            obj.dbn.fit(Data,Label);
            % Freeze learning of DBN and train the linear perceptron.
            Data        = obj.getDataBlockExclude(ExcludeState);
            Label       = obj.getLabelBlockExclude(ExcludeState); 
            nLabel      = size(Label,1); % Convert labels to 0, 1 and 1D.
            [~, Label]  = max(Label,[],2);
            Label       = Label - 1;
            A           = obj.dbn.probe(Data); % nLayer x nHidden x 16
            Wc          = obj.dbn.getLastLayer.getWeightForClass(); % 2 x nHidden
            Data        = sum(repmat(Wc(1,:),[nLabel 1])'.*squeeze(A(end,:,:)))';
            obj.lp.train(Data,Label); % Requires label numbers 0 and 1.
        end
        function err = testError(obj)
            Data        = obj.getDataBlock();
            Label       = obj.getLabelBlock();
            nLabel      = size(Label,1);
            A           = obj.dbn.probe(Data); % nLayer x nHidden x 16
            Wc          = obj.dbn.getLastLayer.getWeightForClass(); % 2 x nHidden
            Data        = sum(repmat(Wc(1,:),[nLabel 1])'.*squeeze(A(end,:,:)))';
            L           = obj.lp.predict(Data);
            [~, Label]  = max(Label,[],2);
            Label       = Label - 1;
            err         = sum(L~=Label)/length(L);
        end
        function [A Wc] = getDBNActivationSortedByWeights(obj)
            Data    = obj.getDataBlock();
            A       = obj.dbn.probe(Data); % nLayer x nHidden x 16
            Wc      = obj.dbn.getLastLayer.getWeightForClass(); % 2 x nHidden
            [~, Index] = sort(abs(Wc(1,:)),2,'descend');
            A       = A(:,Index,:);
            Wc      = Wc(1,Index);
        end
        function D = getLPData(obj)
            Data    = obj.getDataBlock();
            A       = obj.dbn.probe(Data); % nLayer x nHidden x 16
            Wc      = obj.dbn.getLastLayer.getWeightForClass(); % 2 x nHidden
            D       = sum(repmat(Wc(1,:),[nLabel 1])'.*squeeze(A(end,:,:)))';
        end
        function obj = setBlockTrain(obj,flag)
            obj.blockTrain = flag;
            obj.dbn.setBlockTrain(flag);
        end
    end
    methods (Static = true)
        function id = getIdentifier()
            id = 'DBNaLP';
        end
    end
end