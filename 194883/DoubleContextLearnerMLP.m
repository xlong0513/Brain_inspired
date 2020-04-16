classdef DoubleContextLearnerMLP < DoubleContextLearner 
    % DoubleContextLearner
    % Uses a Multi Layer Perceptron (MLP) network to learn the 
    % double-context task.
    %   Florian Raudies, 01/30/2014, Boston University.
    properties (SetAccess = private)
        mlp         % Multi layer perceptron network.
    end
    methods
        function obj = DoubleContextLearnerMLP(LetterLabel,NumberLabel,nHidden)
            obj = obj@DoubleContextLearner(LetterLabel,NumberLabel);
            obj.mlp = MultiLayerPerceptronNetwork(nHidden);
        end
        function obj = learn(obj,nBlock,ExcludeState)
            [Data Label] = obj.generateData(nBlock, ExcludeState);
            obj.mlp.train(Data,Label);
        end
        function err = testError(obj)
            Data    = obj.getDataBlock();
            Label   = obj.getLabelBlock();
            L       = obj.mlp.retrieve(Data);
            err     = sum(sum((L>.5)~=Label))/(2*length(L));
        end
        function obj = setBlockTrain(obj,flag)
            obj.blockTrain = flag;
            obj.mlp.setBlockTrain(flag);
        end
        function E = getTrainingError(obj)
            E = obj.mlp.getTrainingError();
        end
    end
    methods (Static = true)
        function id = getIdentifier()
            id = 'MLP';
        end
    end
end