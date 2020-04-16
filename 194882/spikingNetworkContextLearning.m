function [PerCorrect FiringRate FIndex W12perTrial W23perTrial RasterPlot] = ...
    spikingNetworkContextLearning(nTrial)
% spikingNetworkContextLearning
%   nTrial  - Number of trials. Numerous other parameters are specified in 
%             the function itself.
%
% RETURN
%   PerCorrect  - Percent correct with dimensions: nTrial x 1.
%   FiringRate  - Firing rate with dimensions: nTrial x nStim x nHippo.
%   FIndex      - Functional index as binary matrix. A strong weight 
%                 connection from the 1st to 2nd layer indicates this 
%                 connection being funcational. This index has the 
%                 dimensions: 1 x nHippo.
%   W12PerTrial - Weight matrix from 1st to 2nd layer for all trials. This 
%                 matrix has the dimensions: nTrial x nInput x nHippo.
%   W23PerTrial - Weight matrix from 2nd to 3rd layer for all trials. This
%                 matrix has the dimensions: nTrial x nHippo x nOutput.
%
% DESCRIPTION
%   This is the main network simulation. It inlcude the
%   - initialization of the network
%   - the spiking simulation for each trial with its phase of the trial and 
%     phase of reply.
%   - during replay the synaptic weights between 1st/2nd and 2nd/3rd layer
%     are adapted.
%

%   Florian Raudies, 09/07/2014, Boston University.

IntvlTrial  = [0 4000]; % ms
IntvlReplay = [0 400];  % ms
dt          = 0.5;      % time step
dtWindow    = 10;       % ms
TimeTrial   = ( IntvlTrial(1)  : dt : IntvlTrial(2)  )';
TimeReplay  = ( IntvlReplay(1) : dt : IntvlReplay(2) )';
nTimeTrial  = length(TimeTrial);
nTimeReplay = length(TimeReplay);
nWindow     = dtWindow/dt;
V_PEAK      = 0;            % in Volts, these are 0 mV.
V_TH        = -50*10^-3;    % in Volts, these are -50 mV.
V_RESET     = -70*10^-3;    % in Volts, these are -70 mV.
ETA         = 10^-6;        % Threshold for equal

%            Context and Place   Odor
% Input for: A1 | B1 | A2 | B2 | X | Y
InVec       = [1 0 0 0 1 0];

% Per definition the rewarded stimuli are: A1X  A2X  B1Y  B2Y.
reward      = @(X) (X(1) && X(5)) || (X(3) && X(5)) ...
                || (X(2) && X(6)) || (X(4) && X(6));
stimulusIndex = @(X) 1*X(1)+2*X(2)+3*X(3)+4*X(4)+4*X(6);

% Output action vector: Dig | Move
OutVec  = [0 0];

% Number of neurons to simulate the hippocampus.
nHippo  = 8;
nIn     = length(InVec);
nOut    = length(OutVec);
nStim   = 8;

% Randomly initialize the synaptic coupling strengths (weights).
% Per random some of these have place cell selectivity.
W12 = rand(nIn,nHippo);
W23 = rand(nHippo,nOut);
W22 = ones(nHippo,nHippo)   -eye(nHippo); % Inhibition weights
W33 = ones(nOut,nOut)       -eye(nOut);

% Number of steps in the history.
nHist       = 2;
InVecHst    = zeros(nHist,nIn);
HippoVecHst = zeros(nHist,nHippo);
OutVecHst   = zeros(nHist,nOut);

% Buffer for spikes within the STDP window.
nMaxSpike   = 10;
T1          = TimeBuffer(nMaxSpike,nIn,dtWindow);
T2          = TimeBuffer(nMaxSpike,nHippo,dtWindow);
T3          = TimeBuffer(nMaxSpike,nOut,dtWindow);

% For performance reasons randomize all indices.
Index       = rand(nTrial,2);

% Percent correct detected.
PerCorrect  = zeros(nTrial,1);
nMaxTime    = 1200;
nMaxSample  = 100;
RasterPlot  = ManySlotBuffer(nStim*nHippo,nMaxSample,nMaxTime);

% Set the options for the LIF neuron and STDP rule.
opt.V_PEAK  = V_PEAK;
opt.V_TH    = V_TH;
opt.V_RESET = V_RESET;
opt.dt      = dt;

% Define matrices for weights per trial.
W12perTrial = zeros(nTrial,nIn,nHippo);
W23perTrial = zeros(nTrial,nHippo,nOut);

% *************************************************************************
% Loop over all trials.
% *************************************************************************
for iTrial = 1:nTrial,
    % Start trial in a random state.
    InVec = zeros(1,nIn);
    InVec(1+round(Index(iTrial,1)*3))   = 1;
    InVec(1+4+(Index(iTrial,2)>0.5))    = 1;
    % Reset counter for buffers.
    nHistCount  = 0;
    InVecHst(1+nHistCount,:) = InVec;
    TraceV1     = nan(nTimeTrial,nIn);
    TraceV2     = nan(nTimeTrial,nHippo);
    TraceV3     = nan(nTimeTrial,nOut);
    rewarded    = 0;
    nMoveSpike  = 0;
    nDigSpike   = 0;
    % Initialize membrane potentials.
    V1      = repmat(V_RESET,[1 nIn]);
    V2      = repmat(V_RESET,[1 nHippo]);
    V3      = repmat(V_RESET,[1 nOut]);
    NoiseV1 = 10^-6*randn(nTimeTrial,1);
    NoiseV2 = 10^-6*randn(nTimeTrial,1);
    NoiseV3 = 10^-6*randn(nTimeTrial,1);
    nThDigSpike     = 5;
    nThMoveSpike    = 5;
    iLastTime       = 1;
    for iTime = 1 : nTimeTrial,
        t           = TimeTrial(iTime);
        opt.I       = InVec;
        opt.G       = repmat(.1,[1 nIn]);        
        V1          = lifModel(t, V1,opt) + NoiseV1(iTime);
        [~,mi]      = max((V1-V_RESET)*W12 - (V2-V_RESET)*W22);
        opt.I       = zeros(1,nHippo);
        opt.I(mi)   = .98; % nA
        V2          = lifModel(t, V2,opt) + NoiseV2(iTime);
        [~,mi]      = max((V2-V_RESET)*W23 - (V3-V_RESET)*W33);
        opt.I       = zeros(1,nOut);
        opt.I(mi)   = .96; % nA
        V3          = lifModel(t, V3,opt) + NoiseV3(iTime);
        TraceV1(iTime,:) = V1;
        TraceV2(iTime,:) = V2;
        TraceV3(iTime,:) = V3;
        % Register any spikes at the output in the output vector.
        OutVec = zeros(1,nOut);
        OutVec(abs(V3-V_PEAK)<=ETA) = 1;
        % Keep the history of the states/firings.
        if any(abs(V2-V_PEAK)<=ETA)
            HippoVecHst(1+nHistCount,:) = double(abs(V2-V_PEAK)<=ETA);
        end
        if any(abs(V3-V_PEAK)<=ETA)
            OutVecHst(1+nHistCount,:) = double(abs(V3-V_PEAK)<=ETA);
        end
        nMoveSpike  = nMoveSpike    + OutVec(2);
        nDigSpike   = nDigSpike     + OutVec(1);
        % Dig ?
        if nDigSpike>=nThDigSpike,
            OutVec(1)       = 1;
            OutVecHst(1+nHistCount,:) = OutVec;
            rewarded = reward(InVec);
            for iHippo = 1:nHippo,
                iSlot   = sub2ind([nStim nHippo],stimulusIndex(InVec),iHippo);
                DataRow = [iTrial; iTime-iLastTime+1; ...
                    abs(TraceV2(iLastTime:iTime,iHippo)-V_PEAK)<=ETA];
                RasterPlot.addEntryToSlot(iSlot,DataRow);
            end
            break;
        end
        % Move?
        if nMoveSpike>=nThMoveSpike,
            nThMoveSpike    = 5;
            nThDigSpike     = max(nThDigSpike - 1,  0);
            nMoveSpike      = 0;
            OutVec(2)       = 1;
            OutVecHst(1+nHistCount,:) = OutVec;
            for iHippo = 1:nHippo,
                iSlot   = sub2ind([nStim nHippo],stimulusIndex(InVec),iHippo);
                DataRow = [iTrial; iTime-iLastTime+1; ...
                    abs(TraceV2(iLastTime:iTime,iHippo)-V_PEAK)<=ETA];
                RasterPlot.addEntryToSlot(iSlot,DataRow);
            end
            iLastTime = iTime+1;
            % Move to the other place.
            Tmp = InVec(3:4);
            InVec(3:4) = InVec(1:2);
            InVec(1:2) = Tmp;
            % Then percept changes too.
            InVec(5) = InVec(6);
            InVec(6) = ~InVec(5);
            % Increment the counter for the buffer.
            nHistCount = mod(nHistCount + 1,nHist);
            InVecHst(1+nHistCount,:) = InVec;
            % Assume there is a break and all the membrane potential return
            % to their resting state.
            V1 = repmat(V_RESET,[1 nIn]);
            V2 = repmat(V_RESET,[1 nHippo]);
            V3 = repmat(V_RESET,[1 nOut]);            
        end
    end
    PerCorrect(iTrial) = rewarded;
    
    % Replay the sequence with the last 1+nHistCount steps.
    for iHist = 1 : (1+nHistCount),
        InVec       = InVecHst(iHist,:);
        HippoVec    = HippoVecHst(iHist,:);
        OutVec      = OutVecHst(iHist,:);
        % Assume there was a break and all membrane potentials return
        % to their resting state value.
        V1          = repmat(V_RESET,[1 nIn]);
        V2          = repmat(V_RESET,[1 nHippo]);
        V3          = repmat(V_RESET,[1 nOut]);
        TraceV1     = nan(nTimeReplay,nIn);
        TraceV2     = nan(nTimeReplay,nHippo);
        TraceV3     = nan(nTimeReplay,nOut);
        TraceW12    = nan(nTimeReplay,nIn,nHippo);
        TraceW23    = nan(nTimeReplay,nHippo,nOut);
        % Clear out any remaining spike times.
        T1.clear();
        T2.clear();
        T3.clear();
        % Start the replay sequence.
        for iTime = 1 : nTimeReplay,
            t       = TimeReplay(iTime);
            % Replay in forward direction --- per STDP strenghening
            if rewarded,
                opt.I   = InVec;
                V1      = lifModel(t, V1, opt);
                opt.I   = .98*HippoVec;
                V2      = lifModel(t, V2, opt);
                opt.I   = .96*OutVec;
                V3      = lifModel(t, V3, opt);
            % Replay in inverse direction --- per STDP weakening
            else
                opt.I   = .96*InVec;
                V1      = lifModel(t, V1, opt);
                opt.I   = .98*HippoVec;
                V2      = lifModel(t, V2, opt);
                opt.I   = OutVec;
                V3      = lifModel(t, V3, opt);
            end                    
            % Retire spike times which are too old.
            T1.retire(t);
            T2.retire(t);
            T3.retire(t);
            % Register new spike times.
            T1.addTime(t,abs(V1-V_PEAK)<=ETA);
            T2.addTime(t,abs(V2-V_PEAK)<=ETA);
            T3.addTime(t,abs(V3-V_PEAK)<=ETA);
            % Update the synaptic weights.
            if iTime >= nWindow,
                for iIn = 1:nIn,
                    TimePre     = T1.time(iIn);
                    nPre        = length(TimePre);
                    if nPre==0, continue; end
                    for iHippo = 1:nHippo,
                        TimePost    = T2.time(iHippo);
                        nPost       = length(TimePost);
                        % Are there any spikes for the pre- and the 
                        % post-synaptic neuron in the time window?
                        if nPre>0 && nPost>0,
                            opt.TimePre     = TimePre;
                            opt.TimePost    = TimePost;
                            W12(iIn,iHippo) = stdpModel(...
                                                t,W12(iIn,iHippo),opt);
                        end
                    end
                end
                for iHippo = 1:nHippo,
                    TimePre     = T2.time(iHippo);
                    nPre        = length(TimePre);
                    if nPre==0, continue; end
                    for iOut = 1:nOut,
                        TimePost    = T3.time(iOut);
                        nPost       = length(TimePost);
                        % Are there any spikes for the pre- and the 
                        % post-synaptic neuron in the time window?
                        if nPre>0 && nPost>0,
                            opt.TimePre     = TimePre;
                            opt.TimePost    = TimePost;
                            W23(iHippo,iOut) = stdpModel(...
                                                t,W23(iHippo,iOut),opt);
                        end
                    end
                end
            end
            TraceV1(iTime,:)    = V1;
            TraceV2(iTime,:)    = V2;
            TraceV3(iTime,:)    = V3;
            TraceW12(iTime,:,:) = W12;
            TraceW23(iTime,:,:) = W23;
        end
    end
    W12perTrial(iTrial,:,:) = W12;
    W23perTrial(iTrial,:,:) = W23;
end

opt.nTrial      = nTrial;
opt.nStim       = nStim;
opt.nCell       = nHippo;
opt.nMaxSample  = nMaxSample;

% Calculate the firing rate for each hippocampal cell, trial, and stimulus. 
FiringRate = rasterPlotToFiringRate(RasterPlot, opt);
% Calculate a binary vector indicating whether a hippocampal cell is part
% of the functional network or not.
FIndex = max(W12,[],1)>(1-ETA);

