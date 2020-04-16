classdef LinearPerceptron < handle
    % Linear perceptron learns a (D-1) dimensional hyperplane for
    % D-dimensional data. Can only solve linear separable problems with two
    % classes.
    %   Florian Raudies, 01/09/2014, Boston University.
    properties (SetAccess = private)
        W       % Weight.
        theta   % Threshold.
    end
    properties
        my          % Damping parameter.
        nMaxStep    % Maximum number of steps.
        nStep       % Steps taken to converge.
        eta         % Error threshold.
        verbose     % Print infos.
    end
    methods
        % Constructor.
        function obj = LinearPerceptron()
            obj.my          = 1e-3;
            obj.nMaxStep    = 1e+3;
            obj.eta         = 1e-2;
            obj.verbose     = 0;
        end
        % Training algorithm.
        function obj = train(obj, Data, Label)
            nDim        = size(Data,2);
            nData       = size(Data,1);
            obj.W       = rand(nDim,1);
            obj.theta   = rand;
            obj.nStep   = obj.nMaxStep;
            for iStep = 1:obj.nMaxStep
                X       = Data*obj.W;
                errSum  = 0;
                for iData = 1:nData,
                    y           = double( (X(iData) - obj.theta) >= 0);
                    d           = Label(iData) - y;
                    obj.W       = obj.W + obj.my*Data(iData,:)'*d;
                    obj.theta   = obj.theta - obj.my*d;
                    errSum      = errSum + d^2;
                end
                if errSum < obj.eta,
                    obj.nStep = iStep;
                    if obj.verbose, fprintf('%d steps.\n',iStep); end
                    break; 
                end
            end
        end
        % Prediction 'P' of classes for data 'Data'.
        function P = predict(obj,Data)
            P = double( (Data*obj.W - obj.theta) >= 0);
        end
        % Get weight 'W'.
        function W = getWeight(obj)
            W = obj.W;
        end
        % Get threshold 'theta'.
        function theta = getThreshold(obj)
            theta = obj.theta;
        end
    end
end