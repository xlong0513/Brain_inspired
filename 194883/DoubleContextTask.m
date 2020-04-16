classdef DoubleContextTask < handle
    % DoubleContextTask
    % The double context task requires the assocation of 16 stimulus 
    % (A,B,C,D) - context (1,2,3,4) pairs with one of the two responses X
    % or Y.
    %
    % The task is as follows.
    %   ----------------                 -----------
    %   | A1 B1 | A2 B2 |               | X X | Y Y |
    %   | C1 D1 | C2 D2 |  associate    | Y Y | X X |
    %   ----------------   --------->    -----------
    %   | A3 B3 | A4 B4 |               | Y Y | X X |
    %   | C3 D3 | C4 D4 |               | X X | Y Y |
    %   ----------------                 -----------
    %
    %
    %   Florian Raudies, 01/30/2014, Boston University.
    properties (SetAccess = protected)
        LetterLabel 
        NumberLabel
        StateName
        DataBlock
        LabelBlock
        blockTrain  % Train with ordered blocks.
    end    
    methods
        % For the double-conext task call with
        % LetterLabel = {'A','B','C','D'} and 
        % NumberLabel = {'1','2','3','4'}
        function obj = DoubleContextTask(LetterLabel,NumberLabel)
            obj.LetterLabel = LetterLabel;
            obj.NumberLabel = NumberLabel;
            nLetter = length(obj.LetterLabel);
            nNumber = length(obj.NumberLabel);
            obj.StateName   = cell(nLetter * nNumber, 1);
            obj.DataBlock   = zeros(nLetter*nNumber,nLetter+nNumber);
            LabelIndex  = zeros(nLetter*nNumber,1);
            for iLetter = 1:nLetter,
                letter = obj.LetterLabel{iLetter};
                for iNumber = 1:nNumber,
                    iData = sub2ind([nNumber nLetter],iNumber,iLetter);
                    obj.StateName{iData} = [letter, ...
                                            obj.NumberLabel{iNumber}];
                    if iLetter <= nLetter/2,
                        LabelIndex(iData) = iNumber==2 || iNumber==3;
                    else
                        LabelIndex(iData) = ~(iNumber==2 || iNumber==3);
                    end
                    obj.DataBlock(iData,iLetter) = 1;
                    obj.DataBlock(iData,nLetter+iNumber) = 1;
                end
            end
            LabelIndex      = 1 + double(LabelIndex);
            LabelIndex      = sub2ind([nLetter*nNumber 2],...
                                      (1:nLetter*nNumber)',LabelIndex);
            obj.LabelBlock  = zeros(nLetter*nNumber,2);
            obj.LabelBlock(LabelIndex) = 1;
            obj.blockTrain  = 0;
        end
        function [Data Label] = generateData(obj, nBlock, ExcludeState)
            [~, Exclude] = ismember(ExcludeState, obj.StateName);
            Include = setdiff(1:size(obj.DataBlock,1), Exclude);
            Data    = obj.DataBlock(Include,:);
            Label   = obj.LabelBlock(Include,:);
            if ~obj.blockTrain
                Data    = repmat(Data,nBlock,1);
                Label   = repmat(Label,nBlock,1);
                Index   = randperm(length(Label));
                Data    = Data(Index,:);
                Label   = Label(Index,:);
            else
                Index   = arrangeBlocks(size(Data,1),nBlock,1);
                Data    = Data(Index,:);
                Label   = Label(Index,:);
            end
        end
        function Data = getDataBlock(obj)
            Data    = obj.DataBlock;
        end
        function Label = getLabelBlock(obj)
            Label   = obj.LabelBlock;
        end
        function Data = getDataBlockExclude(obj, ExcludeState)
            [~, Exclude] = ismember(ExcludeState, obj.StateName);
            Include = setdiff(1:size(obj.DataBlock,1), Exclude);
            Data    = obj.DataBlock(Include,:);
        end
        function Label = getLabelBlockExclude(obj, ExcludeState)
            [~, Exclude] = ismember(ExcludeState, obj.StateName);
            Include = setdiff(1:size(obj.DataBlock,1), Exclude);
            Label   = obj.LabelBlock(Include,:);
        end
    end
end