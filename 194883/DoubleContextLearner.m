classdef DoubleContextLearner < DoubleContextTask
    % DoubleContextLearner
    % Abstract class for any class that provides learning of the
    % double-context task.
    %   Florian Raudies, 01/30/2014, Boston University.
    methods
        % Constructor.
        function obj = DoubleContextLearner(LetterLabel,NumberLabel)
            obj = obj@DoubleContextTask(LetterLabel,NumberLabel);
        end
    end
    methods (Abstract = true)
        % Learning of the task with the number of blocks nBlock
        % (repititions of the original data) while excluding any
        % stimulus-context combinations listed in ExcludeState.
        obj = learn(obj,nBlock,ExcludeState)
        % Test all stimulus-context combinations and calculate the error
        % rate.
        err = testError(obj)
        % Train the stimulus-context combinations in the same order rather
        % than random order in the different epochs.
        obj = setBlockTrain(obj,flag)
    end
    methods (Abstract = true, Static = true)
        % Get a unique identifier for the learner.
        id = getIdentifier()
    end
end