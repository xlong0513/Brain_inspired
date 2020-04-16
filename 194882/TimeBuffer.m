classdef TimeBuffer < handle
    % Time buffer -- Is a single buffer for time stamps. The method retire
    % removes all entries which are older than t-nHistory.
    
    %   Florian Raudies, 09/07/2014, Boston University.
    properties
        nEntry      % Number of entries.
        nBuffer     % Number of buffers.
        nHistory    % Time in steps to go back.
        Counter     % Counter for each buffer.
        Buffer      % These are the buffers that hold the time values.
    end
    methods
        % Constructor.
        function obj = TimeBuffer(nEntry,nBuffer,nHistory)
            obj.nEntry      = nEntry;
            obj.nBuffer     = nBuffer;
            obj.nHistory    = nHistory;
            obj.Counter     = zeros(nBuffer,1);
            obj.Buffer      = zeros(nEntry,nBuffer);
        end
        function obj = clear(obj)
            obj.Counter = zeros(obj.nBuffer,1);
            obj.Buffer  = zeros(obj.nEntry,obj.nBuffer);
        end
        % Retire entries.
        function obj = retire(obj,time)
            for iBuffer = 1:obj.nBuffer,
                if obj.Counter(iBuffer)==0, continue; end
                ToRetire = obj.Buffer(1:obj.Counter(iBuffer),iBuffer) < time-obj.nHistory;
                if sum(ToRetire)==0, continue; end
                Index = sum(ToRetire)+1 : obj.Counter(iBuffer);
                if isempty(Index), obj.Counter(iBuffer) = 0; continue; end
                obj.Buffer(1:length(Index),iBuffer) = obj.Buffer(Index,iBuffer);
                obj.Counter(iBuffer) = length(Index);
            end
        end
        % Add time.
        function obj = addTime(obj,time,ToBuffer)
            if any(ToBuffer)
                for iToBuffer = find(ToBuffer),
                    obj.Buffer(1+obj.Counter(iToBuffer),iToBuffer) = time;
                    obj.Counter(iToBuffer) = obj.Counter(iToBuffer) + 1;
                end
            end
        end
        % Retrieve time for iBuffer.
        function Time = time(obj,iBuffer)
            Time = obj.Buffer(1:obj.Counter(iBuffer),iBuffer);
        end
        % Print the contents of the buffer.
        function print(obj)
            for iBuffer = 1:obj.nBuffer,
                fprintf('buffer %d: ',iBuffer);
                for iEntry = 1:obj.Counter(iBuffer),
                    fprintf('%d, ',obj.Buffer(iEntry,iBuffer));
                end
                fprintf('\n');
            end            
        end
    end
end