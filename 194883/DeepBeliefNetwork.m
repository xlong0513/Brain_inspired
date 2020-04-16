classdef DeepBeliefNetwork < handle
    % A deep belief network composed of 'nLayer' layers of restricted
    % Boltzmann machines.
    %
    %   Florian Raudies, 01/08/2014, Boston University.
    %
    %   This is a re-implementation of Andrej Karpathy's code which was 
    %   based on Kevin Swersky and Ruslan Salakhutdinov's code.
    properties (SetAccess = private)
        Layer        % Layers for Restricted Boltzmann Machines.
    end
    properties
        DimLayer     % Dimensions for each layer.
        nLayer       % Number of layers.
    end
    methods
        % Constructor.
        function obj = DeepBeliefNetwork(DimLayer)
            obj.DimLayer = DimLayer;
            obj.nLayer   = length(DimLayer);
            obj.Layer    = cell(obj.nLayer,1);
            for iLayer = 1:obj.nLayer,
                obj.Layer{iLayer} = RestrictedBoltzmannMachine();
            end
        end
        % Fit 'Data' to 'Labels' through a deep belief network.
        function obj = fit(obj, Data, Label)
            if obj.nLayer > 1,
                obj.Layer{1}.trainContrastiveConvergence(...
                    Data, obj.DimLayer(1));
                for iLayer = 2:obj.nLayer-1
                    obj.Layer{iLayer}.trainContrastiveConvergence(...
                        obj.Layer{iLayer-1}.T, obj.DimLayer(iLayer));
                end
                % Train the last layer with labels.
                obj.Layer{obj.nLayer}.fitContrastiveConvergence(...
                        obj.Layer{obj.nLayer-1}.T, Label, ...
                        obj.DimLayer(obj.nLayer));
            else
                obj.Layer{1}.fitContrastiveConvergence(...
                        Data, Label, obj.DimLayer(1));
            end
        end
        % Predict labels depending on 'Data'.
        function P = predict(obj, Data)
            for iLayer = 1:obj.nLayer-1,
                Data = obj.Layer{iLayer}.visibleToHidden(Data);
            end
            P = obj.Layer{obj.nLayer}.predict(Data);
        end
        % Retrieve the internal activation when probing with 'Data'.
        % This method assumes that all hidden layers have the same size.
        function P = probe(obj, Data)
            nState = size(Data,1);
            P = zeros(obj.nLayer,obj.DimLayer(1),nState);
            for iState = 1:nState,
                Data4State = Data(iState,:);
                for iLayer = 1:obj.nLayer,
                    Data4State = obj.Layer{iLayer}.visibleToHidden(Data4State);
                    P(iLayer,:,iState) = Data4State;
                end
            end
        end
        % Get network representation for a 'iLayer'.
        function rbm = getLayer(obj,iLayer)
            rbm = obj.Layer{iLayer};
        end
        % Get network representation for the last layer.
        function rbm = getLastLayer(obj)
            rbm = obj.Layer{obj.nLayer};
        end
        % Set flag which indicates training in ordered blocks.
        function obj = setBlockTrain(obj,flag)
            for iLayer = 1:obj.nLayer, 
                obj.Layer{iLayer}.setBlockTrain(flag);
            end
        end
    end
end