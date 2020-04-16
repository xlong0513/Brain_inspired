classdef MultiLayerPerceptronNetwork < handle
    % Multi layer perceptron network (here two layers).
    %
    %   Florian Raudies, 01/30/2014, Boston University.
    properties (SetAccess = private)
        W1          % Weights between input and hidden layer.
        W2          % Weights between hidden and output layer.
        Theta1      % Thresholds for hidden layer.
        Theta2      % Thresholds for output layer.
        E           % Training error.
        blockTrain  % No shuffeling of data when training.
    end
    properties
        nDimIn
        nDimOut
        nHidden     % Number of hidden neurons.
        nMaxEpoch   % Maximum epoch number.
        beta        % Parameter of Fermi function.
        eta         % Learning rate.
        fSgd        % Sigmoid function.
        fSgdDiff    % Derivative of sigmoid function.
    end
    methods
        % Constructor sets default values and initializes network weights
        % and thresholds to allow for incremental or repeated learning.
        function obj = MultiLayerPerceptronNetwork(nHidden)
            obj.nHidden     = nHidden;
            obj.nMaxEpoch   = 50;
            obj.beta        = 1;
            obj.eta         = 0.1;
            obj.fSgd        = @(X) 1./(1 + exp(-obj.beta*X));
            obj.fSgdDiff    = @(X) obj.beta*obj.fSgd(X).*(1-obj.fSgd(X));
        end
        function obj = train(obj,Data,Label)
            nData       = size(Data,1);
            obj.nDimIn  = size(Data,2);
            obj.nDimOut = size(Label,2);
            obj.W1      = rand(obj.nDimIn,obj.nHidden)   - 0.5;
            obj.W2      = rand(obj.nHidden,obj.nDimOut)  - 0.5;
            obj.Theta1  = zeros(obj.nHidden,1);
            obj.Theta2  = zeros(obj.nDimOut,1);
            if ~obj.blockTrain
                Shuffle     = randperm(nData);
                Data        = Data(Shuffle,:);
                Label       = Label(Shuffle,:);
            end
            obj.E       = zeros(obj.nMaxEpoch,1);
            for iEpoch = 1:obj.nMaxEpoch,
                msError = 0;
                for iData = 1:nData,
                    X = Data(iData,:)';
                    L = Label(iData,:)';
                    % Forward sweep through the MLP.
                    U1 = obj.W1'*X - obj.Theta1;
                    Y1 = obj.fSgd(U1);
                    U2 = obj.W2'*Y1 - obj.Theta2;
                    Y2 = obj.fSgd(U2);
                    % Calculate the output error.
                    D2 = (L-Y2).*obj.fSgdDiff(U2);
                    % Propagate the output error backward through the
                    % network.
                    D1 = (obj.W2*D2).*obj.fSgdDiff(U1);
                    % Learn by adaptation of weights and thresholds.
                    obj.W2      = obj.W2 + obj.eta*Y1*D2';
                    obj.W1      = obj.W1 + obj.eta*X*D1';
                    obj.Theta2  = obj.Theta2 - obj.eta*D2;
                    obj.Theta1  = obj.Theta1 - obj.eta*D1;
                    msError     = msError + sum((L-Y2).^2);
                end
                obj.E(iEpoch) = sqrt(msError/nData);
            end
        end
        function obj = setBlockTrain(obj,flag)
            obj.blockTrain = flag;
        end
        function L = retrieve(obj,Data)
            nData = size(Data,1);
            L = zeros(nData,obj.nDimOut);
            for iData = 1:nData,
                X = Data(iData,:)';
                % Forward sweep through the MLP.
                U1 = obj.W1'*X - obj.Theta1;
                Y1 = obj.fSgd(U1);
                U2 = obj.W2'*Y1 - obj.Theta2;
                L(iData,:) = obj.fSgd(U2)';
            end            
        end
        function E = getTrainingError(obj)
            E = obj.E;
        end
    end
end