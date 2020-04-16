classdef DoubleContextLearnerDBN < DoubleContextLearner 
    % DoubleContextLearner
    % Uses a Deep Belief Network (DBN) to learn the double-context task.
    %   Florian Raudies, 01/30/2014, Boston University.
    properties (SetAccess = private)
        dbn         % Deep belief network.
    end
    methods
        function obj = DoubleContextLearnerDBN(...
                LetterLabel,NumberLabel,nHidden,nLayer)
            obj = obj@DoubleContextLearner(LetterLabel,NumberLabel);
            obj.dbn = DeepBeliefNetwork(repmat(nHidden,[1 nLayer]));
        end
        function obj = learn(obj,nBlock,ExcludeState)
            [Data Label] = obj.generateData(nBlock, ExcludeState);
            obj.dbn.fit(Data,Label);
        end
        function err = testError(obj)
            Data    = obj.getDataBlock();
            Label   = obj.getLabelBlock();
            L       = obj.dbn.predict(Data);
            err     = sum(sum(L~=Label))/(2*length(L));
        end
        function obj = setBlockTrain(obj,flag)
            obj.blockTrain = flag;
            obj.dbn.setBlockTrain(flag);
        end
    end
    methods (Static = true)
        function id = getIdentifier()
            id = 'DBN';
        end
    end
end