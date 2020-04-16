function[in,pos,LFP,gridCells,cycle,cycleDecoding] = humanGridSim(in)
%% Script to simulate grid cell activity in a 1D or 2D environment using
%  a phenomenological model adapted from Chadwick et al. eLife 4:e03542 
%  (2015), and then decode movement trajectory from grid cell firing rate 
%  and phase in each oscillatory cycle. 
%  Daniel Bush, UCL (2019) drdanielbush@gmail.com
%
%  Described in Bush and Burgess (2020) Detection and advantages of phase
%  coding in the absence of rhythmicity. Hippocampus (in press)
%  https://doi.org/10.1002/hipo.23199
%
%  NOTE:
%  Requires Tom O'Haver's fast smoothing function from:
%  https://www.mathworks.com/matlabcentral/fileexchange/19998-fast-smoothing-function
%
%  INPUTS, as fields of the 'in' structure, default values indicated in []:
%  gridScales   = Grid scales (cm), for five modules [30*1.4.^(0:4)]
%  nGCs      	= Number of grid cells per module [40]
%  sampleRate   = Simulation time step (Hz) [200]
%  sigmas       = Firing field width parameter (cm) [in.gridScales./10]
%  meanRate     = Mean firing rate (Hz) [1]
%  phaseLock    = Phase locking (true) or precession (false)? [false]
%  phaseMod     = Extent of phase coding 'k' (au) [1.5]
%  speedSlope   = Slope of running speed v firing rate (Hz/cm/s) [5/30]
%  environment  = Environment type (1D or 2D) ['2D']
%  trackFile    = 2D tracking data file [ceil(rand*3)]
%  vRange       = Moving speed range for 1D environments [2 30]
%  lfpType   	= LFP type ('constant' or 'human') ['human']
%  freqRange   	= LFP frequency range for human LFP (Hz) [2 20]
%  nPhaseBins   = Number of phase bins for decoding analysis [5]
%  fieldRateVar	= Standard deviation of in-field firing variability [0]



%% Provide some general parameters for the analysis
in              = dealWithInputs(in);

%  ...generate or load some tracking data
disp('Generating / importing tracking data...'); drawnow
pos             = trackingData(in);

%  ...generate or load some LFP data
disp('Generating / importing LFP data...'); drawnow
LFP             = lfpData(in,pos);

%  ...simulate grid cell firing patterns
disp('Generating grid cell spike trains...'); drawnow
gridCells       = gridCellFiring(in,pos,LFP);

%  ...chunk the data by oscillatory cycle
disp('Chunking the data by oscillatory cycle...'); drawnow
cycle           = chunkData(in,pos,LFP,gridCells);

%  ...decode location / movement direction / running speed / anxiety
disp('Decoding movement trajectory from activity in each cycle...'); drawnow
cycleDecoding	= decodeCycles(in,pos,gridCells,cycle);

%  ...remove superfluous variables to save memory
gridCells       = rmfield(gridCells,'r_rate');
gridCells       = rmfield(gridCells,'r_phase');
if isfield(gridCells,'r_pure')
    gridCells  	= rmfield(gridCells,'r_pure');
end

%  ...and plot some summary decoding figures
plotResults(cycle,cycleDecoding);
clc

end



%% Function to organise input settings
function[in]	= dealWithInputs(in)
if ~isfield(in,'gridScales') || isempty(in.gridScales)
    in.gridScales   = 30*1.4.^(0:4);        % Grid scales (cm), five modules
end
if ~isfield(in,'nGCs') || isempty(in.nGCs)
    in.nGCs      	= 40;                   % Number of grid cells per module
end
if ~isfield(in,'sampleRate') || isempty(in.sampleRate)
    in.sampleRate   = 200;                  % Simulation time step (Hz)
end
if ~isfield(in,'sigmas') || isempty(in.sigmas)
    in.sigmas       = in.gridScales./10;    % Firing field width parameter (cm)
end
if ~isfield(in,'meanRate') || isempty(in.meanRate)
    in.meanRate     = 1;                    % Mean firing rate (Hz)
end
if ~isfield(in,'phaseLock') || isempty(in.phaseLock)
    in.phaseLock    = false;                % Phase locking (true) or precession (false)?
end
if ~isfield(in,'phaseMod') || isempty(in.phaseMod)
    in.phaseMod     = 1.5;                  % Extent of phase coding 'k' (au)
end
if ~isfield(in,'speedSlope') || isempty(in.speedSlope)
    in.speedSlope   = 5/30;                 % Slope of running speed v firing rate (Hz/cm/s)
end
if ~isfield(in,'environment') || isempty(in.environment)
    in.environment  = '2D';                 % Environment type (1D or 2D)
end
if strcmp(in.environment,'2D') && (~isfield(in,'trackFile') || isempty(in.trackFile))
    in.trackFile    = ceil(rand*3);         % 2D tracking data file
end
if ~isfield(in,'vRange') || isempty(in.vRange)
    in.vRange       = [2 30];               % Moving speed range for 1D environments
end
if ~isfield(in,'lfpType') || isempty(in.lfpType)
    in.lfpType   	= 'human';              % LFP type ('constant' or 'human')
end
if ~isfield(in,'freqRange') || isempty(in.freqRange)
    in.freqRange   	= [2 20];               % LFP frequency range for human LFP (Hz)
end
if ~isfield(in,'nPhaseBins') || isempty(in.nPhaseBins)
    in.nPhaseBins   = 5;                    % Number of phase bins for decoding analysis
end
if ~isfield(in,'fieldRateVar') || isempty(in.fieldRateVar)
    in.fieldRateVar = 0;                    % Standard deviation of in-field firing variability
end
end



%% Function to generate or load tracking data
function[pos]   = trackingData(in)
switch in.environment
    
    % ...either randomly generate 1D tracking data
    case '1D'        
        t_log   = 1/in.sampleRate : 1/in.sampleRate : 30;   % Create a time base (s)
        v_log   = cumsum(randn(length(t_log),1));           % Assign time varying velocity (cm/s)
        v_log   = ((v_log-min(v_log))./range(v_log).*diff(in.vRange)+in.vRange(1))'; clear v_range
        x_log   = cumsum(v_log ./ in.sampleRate);           % Compute x co-ordinates
        y_log   = ones(size(x_log));                        % Assign y co-ordinates        
        
	% ...randomly generate longer 1D tracking data
    case '1Dlong'        
        t_log   = 1/in.sampleRate : 1/in.sampleRate : 300; 	% Create a time base (s)
        v_log   = cumsum(randn(length(t_log),1));           % Assign time varying velocity (cm/s)
        v_log   = ((v_log-min(v_log))./range(v_log).*diff(in.vRange)+in.vRange(1))'; clear v_range
        x_log   = cumsum(v_log ./ in.sampleRate);           % Compute x co-ordinates
        y_log   = ones(size(x_log));                        % Assign y co-ordinates        
        
  	% or load and manipulate 2D tracking data
    case '2D'
        load('trackingData.mat')
        string  = ['xy = Square' int2str(in.trackFile) './PixPerM*100;']; eval(string); clear string Square1 Square2 Square3 PixPerM
        t_log   = linspace(0,size(xy,1)/Fs,size(xy,1)*in.sampleRate/Fs);
        t_dat   = linspace(0,size(xy,1)/Fs,size(xy,1));
        x_log   = xy(:,1);
        y_log   = xy(:,2); clear xy
        v_log   = nan(size(x_log));
        v_log(1:end-1)	= sqrt(diff(x_log).^2 + diff(y_log).^2) * Fs; clear Fs
        v_log(end)      = v_log(end-1);
        x_log   = interp1(t_dat,x_log,t_log);
        y_log   = interp1(t_dat,y_log,t_log);
        v_log   = interp1(t_dat,v_log,t_log); clear t_dat
        
    case '2Dlin'        
        boxSize     = [100 100];    % Box size, cm
        runSpeed    = 30;           % Constant running speed, cm/s
        x_log       = nan(1,ceil(prod(boxSize)/runSpeed*in.sampleRate));
        y_log       = nan(1,ceil(prod(boxSize)/runSpeed*in.sampleRate));
        v_log       = runSpeed*ones(1,ceil(prod(boxSize)/runSpeed*in.sampleRate));
        t_log       = linspace(1/in.sampleRate,ceil(prod(boxSize)/runSpeed),ceil(prod(boxSize)/runSpeed*in.sampleRate));
        x_log(1)    = 1;
        y_log(1)    = 1;
        dir         = 1;
        for t       = 2 : length(x_log)            
            x_log(1,t)  = x_log(1,t-1) + dir.*runSpeed/in.sampleRate;
            y_log(1,t)  = y_log(1,t-1);
            if x_log(1,t) > boxSize(1) || x_log(1,t) < 0
                dir         = dir * -1;
                y_log(1,t)  = y_log(1,t) + 1;
            end            
        end
        clear t dir
        
end
pos.x_log       = x_log'; clear x_log
pos.y_log       = y_log'; clear y_log
pos.v_log       = v_log'; clear v_log
pos.t_log       = t_log'; clear t_log
pos.dt          = 1 / in.sampleRate;

end



%% Function to generate or load LFP data
function[LFP]   = lfpData(in,pos)
switch in.lfpType
    case 'human'                
        load('sampleEEG.mat');
        [b,a]           = butter(2,in.freqRange/(Fs/2));                % Generate second order Butterworth filter
        freq            = filtfilt(b,a,eeg(:,1)); clear b a             % Filter EEG data in frequency range of interest
        freq            = angle(hilbert(freq));                         % Get the phase at each time point
        freq            = angle(exp(1i.*diff(freq)));                   % Get the phase difference and wrap around
        freq(freq<0)    = (in.freqRange(1)*2*pi)/Fs;                    % Ignore negative frequencies
        freq            = fastsmooth(freq,Fs/20,3,1);                   % Smooth the data
        freq            = interp1(eeg(2:end,2), freq, pos.dt : pos.dt : eeg(end,2));
        eeg             = interp1(eeg(:,2), eeg(:,1), pos.dt : pos.dt : eeg(end,2));
        freq            = (freq./(2*pi)).*Fs; clear Fs                  % Convert to dynamic frequency
        startInd        = ceil(rand*(length(freq)-length(pos.t_log))); 	% Choose a random segment
        freq            = freq(startInd:startInd+length(pos.t_log)-1);  % Crop the frequency data
        eeg             = eeg(startInd:startInd+length(pos.t_log)-1);   % Crop the eeg data
        phase           = cumsum([0 2.*pi.*pos.dt.*freq(2:end)]);      	% Compute LFP phase at each timepoint        
        anxiety        	= (freq-min(freq))./range(freq); clear startInd	% Ascribe dynamic 'anxiety' variable
    
    case 'constant'
        anxiety         = ones(size(pos.t_log))';
        freq            = 8 .* anxiety;
        phase           = cumsum([0 2.*pi.*pos.dt.*freq(2:end)]);
        eeg             = cos(phase);
        
end
LFP.eeg         = eeg;
LFP.phase       = phase;
LFP.freq        = freq;
LFP.anxiety     = anxiety;
end



%% Function to generate grid cell firing patterns
function[gridCells]     = gridCellFiring(in,pos,LFP)
nCells                  = length(in.gridScales)*in.nGCs;
switch in.environment
    case {'1D','1Dlong'}
        r_rate          = nan(nCells,length(pos.x_log));   % Assign memory for the rate code of each grid cell
        r_phase         = nan(nCells,length(pos.x_log));   % Assign memory for the phase code of each grid cell
        r_cell          = nan(nCells,length(pos.x_log));   % Assign memory for the combined firing output of each cell
        for module      = 1 : length(in.gridScales)
            for cell    = 1 : in.nGCs
                ind             = (module-1)*in.nGCs + cell;
                x_centres       = ((cell-1)/in.nGCs)*in.gridScales(module)-in.gridScales(module) : in.gridScales(module) : max(pos.x_log) + in.gridScales(module);
                offset          = min(abs(repmat(pos.x_log,1,length(x_centres)) - repmat(x_centres,length(pos.x_log),1)),[],2);
                r_rate(ind,:)   = exp(-(offset.^2) ./ (2.*in.sigmas(module)).^2);
                if in.fieldRateVar > 0
                    if ind      == 1
                       r_pure   = nan(nCells,length(pos.x_log));
                    end
                    r_var       = zeros(size(r_rate(ind,:)));
                    [~,pks]     = findpeaks(offset);
                    r_var(pks)  = 1;
                    r_var       = 1 + cumsum(r_var);
                    fieldRates  = 1+randn(length(pks)+1,1) .* in.fieldRateVar;
                    fieldRates(fieldRates<0)    = 0; clear pks
                    r_var       = fieldRates(r_var)'; clear fieldRates
                    r_pure(ind,:)   = r_rate(ind,:);
                    r_rate(ind,:)   = r_rate(ind,:) .* r_var; clear r_var
                end
                clear offset
                if in.phaseLock
                    offset      = pi*ones(size(pos.x_log));
                else
                    offset     	= pos.x_log - x_centres(1) - in.gridScales(module)/2; clear x_centres
                    offset      = mod((-1/in.gridScales(module)).*offset,1)*2*pi;                    
                end
                r_phase(ind,:) 	= exp(in.phaseMod*cos(offset-LFP.phase')); clear offset
                r_phase(ind,:)  = r_phase(ind,:) ./ max(r_phase(ind,:));
                r_cell(ind,:)   = r_rate(ind,:) .* r_phase(ind,:);
                r_cell(ind,:)   = r_cell(ind,:) .* LFP.freq .* in.speedSlope .* pos.v_log';
                r_cell(ind,:)   = r_cell(ind,:) ./ sum(r_cell(ind,:)) .* in.meanRate .* range(pos.t_log); clear ind
            end
        end
        firingRates     = poissrnd(r_cell);
        clear module cell r_cell sigmas
        
    case {'2D','2Dlin'}
        r_rate          = nan(nCells,length(pos.x_log));   % Assign memory for the rate code of each grid cell
        r_phase         = nan(nCells,length(pos.x_log));   % Assign memory for the phase code of each grid cell
        r_cell          = nan(nCells,length(pos.x_log));   % Assign memory for the combined firing output of each cell
        for module      = 1 : length(in.gridScales)
            for cell    = 1 : in.nGCs
                ind             = (module-1)*in.nGCs + cell;
                xFields         = ceil(range(pos.x_log)/in.gridScales(module)) + 4;
                yFields         = ceil(range(pos.y_log)/(in.gridScales(module)*sind(60))) + 4;
                [x_centres,y_centres] = gridFields2D(cell,in.nGCs,in.gridScales(module),xFields,yFields); clear xFields yFields
                
                centres         = [repmat(reshape(x_centres,[1 1 length(x_centres)]),length(pos.x_log),1) repmat(reshape(y_centres,[1 1 length(y_centres)]),length(pos.y_log),1)];
                [offset,c]      = min(sqrt(sum((repmat([pos.x_log pos.y_log],[1 1 length(x_centres)])-centres).^2,2)),[],3); clear centres
                r_rate(ind,:)   = exp(-(offset.^2) ./ (2.*in.sigmas(module)).^2); clear offset
                if in.phaseLock
                    offset          = pi*ones(size(pos.x_log));
                else
                    offset          = nan(size(pos.x_log));
                    offset(1:end-1)	= distToFieldCentre([pos.x_log pos.y_log],[x_centres(c) y_centres(c)]) - in.gridScales(module)/2;
                    offset          = mod((1/in.gridScales(module)).*offset,1)*2*pi;
                end
                r_phase(ind,:) 	= exp(in.phaseMod*cos(offset-LFP.phase')); clear offset
                r_phase(ind,:)  = r_phase(ind,:) ./ max(r_phase(ind,:));
                if in.fieldRateVar
                    if ind      == 1
                        r_pure	= nan(nCells,length(pos.x_log));
                    end                    
                    fieldRates  = 1+randn(length(x_centres),1) .* in.fieldRateVar;
                    fieldRates(fieldRates<0)    = 0;
                    c           = fieldRates(c)';
                    r_pure(ind,:)   = r_rate(ind,:);
                    r_rate(ind,:)   = r_rate(ind,:) .* c; 
                end
                clear c                
                r_cell(ind,:)   = r_rate(ind,:) .* r_phase(ind,:);
                r_cell(ind,:)   = r_cell(ind,:) .* LFP.freq .* in.speedSlope .* pos.v_log';
                r_cell(ind,:)   = r_cell(ind,:) ./ nansum(r_cell(ind,:)) .* in.meanRate .* range(pos.t_log); clear ind
            end
        end
        firingRates     = poissrnd(r_cell);
        clear module cell r_cell sigmas
        
end
gridCells.r_rate        = r_rate; clear r_rate
gridCells.r_phase       = r_phase; clear r_phase
gridCells.spikeTrains   = firingRates; clear firingRates
if in.fieldRateVar > 0
    gridCells.r_pure    = r_pure; clear r_pure
end
gridCells.spikeTrains(isnan(gridCells.spikeTrains)) = 0;
end



%% Function to compute field centres for arbitrary 2D grid
function[x_centres,y_centres] = gridFields2D(cell,nGCs,scale,xFields,yFields)

%  Assign some memory
x_centres       = nan(xFields*yFields,1);
y_centres       = nan(xFields*yFields,1);

%  Generate a template to work with
xTemplate       = 0 : xFields-1;
yTemplate       = zeros(1,xFields);

%  Then generate a unit grid field of the requisite size
for yShift      = 0 : yFields-1
    inds        = (1:xFields) + yShift*xFields;
    if mod(yShift,2) == 0
        x_centres(inds)     = xTemplate;
    else
        x_centres(inds)     = xTemplate + cosd(60);
    end
    y_centres(inds)         = yTemplate + yShift*sind(60); clear inds
end
clear yShift

%  Rescale
x_centres       = x_centres * scale;
y_centres       = y_centres * scale;

%  Compute the cell dependent offset
xcoord          = mod(cell-1,sqrt(nGCs));
ycoord          = floor((cell-1)./sqrt(nGCs));
x_shift         = (xcoord + ycoord * cosd(60))/ sqrt(nGCs); clear xcoord
y_shift         = ycoord * sind(60) / sqrt(nGCs); clear ycoord

%  Shift to the origin of choice
x_centres       = x_centres - 2*scale + x_shift*scale;
y_centres       = y_centres - 2*scale + y_shift*scale;
end



%% Function to compute linear distance to perpendicular line through field
%  centre along current trajectory
function[dist]  = distToFieldCentre(loc,grid)

heading    	= diff(loc(:,2)) ./ diff(loc(:,1));
intercept   = loc(1:end-1,2) - heading .* loc(1:end-1,1);
perp        = -1./heading;
perpInt     = grid(1:end-1,2) - perp .* grid(1:end-1,1);

crossX      = (intercept - perpInt) ./ (perp - heading); clear perp perpInt
crossY      = heading .* crossX + intercept; clear heading intercept
crossX(isnan(crossX))   = grid(isnan(crossX),1);
crossY(isnan(crossY))   = loc(isnan(crossY),2);

dist        = sqrt((loc(1:end-1,1) - crossX).^2 + (loc(1:end-1,2) - crossY).^2); 
dist2       = sqrt((loc(2:end,1) - crossX).^2   + (loc(2:end,2) - crossY).^2); clear crossX crossY
dist(dist2>dist)    = -dist(dist2>dist);

end



%% Function to split the data into expected and actual firing rate per 
%  oscillatory cycle, and record mean position in each cycle
function[cycle] = chunkData(in,pos,LFP,gridCells)
[~,peaks]                   = findpeaks(cos(LFP.phase));
nCells                      = length(in.gridScales)*in.nGCs;
nCycles                     = length(peaks)-1;
cycle.nSpikes               = nan(nCells,nCycles);
cycle.expRate               = nan(nCells,length(peaks)-1);
cycle.actLoc                = nan(nCycles,1);
cycle.meanSpeed             = nan(nCycles,1);
cycle.meanAnxiety           = nan(nCycles,1);
cycle.cycleLength           = nan(nCycles,1);
cycle.nSpikes_binned        = nan(nCells,in.nPhaseBins,nCycles);
cycle.expRate_binned        = nan(nCells,in.nPhaseBins,nCycles);
cycle.expPhase_binned       = nan(nCells,in.nPhaseBins,nCycles);
cycle.expComb_binned        = nan(nCells,in.nPhaseBins,nCycles);
cycle.actLoc_binned         = nan(in.nPhaseBins,2,nCycles);
cycle.realVec               = nan(nCycles,1);
if in.fieldRateVar>0
    cycle.expPure           = nan(nCells,length(peaks)-1);
    cycle.expPure_binned	= nan(nCells,in.nPhaseBins,nCycles);
end
for c                       = 1 : nCycles
        
    cycle.cycleLength(c,1)  = (peaks(c+1)-peaks(c)).*pos.dt;
    cycle.nSpikes(:,c)      = sum(gridCells.spikeTrains(:,peaks(c):peaks(c+1)),2);
    cycle.expRate(:,c)      = mean(gridCells.r_rate(:,peaks(c):peaks(c+1)),2);
    cycle.actLoc(c,1)       = mean(pos.x_log(peaks(c):peaks(c+1)));
    cycle.actLoc(c,2)       = mean(pos.y_log(peaks(c):peaks(c+1)));
    cycle.meanSpeed(c,1)    = mean(pos.v_log(peaks(c):peaks(c+1)));       
    cycle.meanAnxiety(c,1)  = mean(LFP.anxiety(peaks(c):peaks(c+1)));           
    
    runningSpikes           = cumsum(sum(gridCells.spikeTrains(:,peaks(c):peaks(c+1)),1));
    [~,phaseBin]            = histc(runningSpikes,linspace(0,runningSpikes(end),in.nPhaseBins+1)); clear runningSpikes    
    phaseBin                = phaseBin';
    phaseBin(phaseBin>in.nPhaseBins) = in.nPhaseBins;
    spikes                  = gridCells.spikeTrains(:,peaks(c):peaks(c+1))';
    [x,y]                   = ndgrid(phaseBin,1:size(spikes,2));    
    cycle.nSpikes_binned(:,:,c)     = accumarray([x(:) y(:)],spikes(:))'; clear spikes    
    rate                    = gridCells.r_rate(:,peaks(c):peaks(c+1))';
    cycle.expRate_binned(:,:,c)     = accumarray([x(:) y(:)],rate(:),[in.nPhaseBins in.nGCs*length(in.gridScales)],@mean)'; clear rate
    phase                   = gridCells.r_phase(:,peaks(c):peaks(c+1))';
    cycle.expPhase_binned(:,:,c)    = accumarray([x(:) y(:)],phase(:),[in.nPhaseBins in.nGCs*length(in.gridScales)],@mean)'; clear phase
    comb                    = (gridCells.r_rate(:,peaks(c):peaks(c+1)) .* gridCells.r_phase(:,peaks(c):peaks(c+1)))';
    cycle.expComb_binned(:,:,c)     = accumarray([x(:) y(:)],comb(:),[in.nPhaseBins in.nGCs*length(in.gridScales)],@mean)'; clear comb
    cycle.actLoc_binned(:,1,c)      = accumarray(phaseBin,pos.x_log(peaks(c):peaks(c+1)),[in.nPhaseBins 1],@mean);
    cycle.actLoc_binned(:,2,c)      = accumarray(phaseBin,pos.y_log(peaks(c):peaks(c+1)),[in.nPhaseBins 1],@mean); 
    if in.fieldRateVar>0
        cycle.expPure(:,c) 	= mean(gridCells.r_pure(:,peaks(c):peaks(c+1)),2);
        pure                = gridCells.r_pure(:,peaks(c):peaks(c+1))';
        cycle.expPure_binned(:,:,c)	= accumarray([x(:) y(:)],pure(:),[in.nPhaseBins in.nGCs*length(in.gridScales)],@mean)'; clear pure
        
        comb                = (gridCells.r_pure(:,peaks(c):peaks(c+1)) .* gridCells.r_phase(:,peaks(c):peaks(c+1)))';
        cycle.expPurC_binned(:,:,c)	= accumarray([x(:) y(:)],comb(:),[in.nPhaseBins in.nGCs*length(in.gridScales)],@mean)'; clear comb
    end    
    clear x y phaseBin
    
    realTraj        = cycle.actLoc_binned(:,:,c);
    bx              = regress(realTraj(:,1),[1:in.nPhaseBins ; ones(1,in.nPhaseBins)]');
    by              = regress(realTraj(:,2),[1:in.nPhaseBins ; ones(1,in.nPhaseBins)]');
    cycle.realVec(c,1)              = atan2(by(1),bx(1)); clear bx by realTraj
end
clear c

%  Normalise each binned rate function by the mean firing rate of that cell
norm             	= sum(cycle.nSpikes,2) ./ sum(cycle.expRate,2); norm(isnan(norm)) = 0;
cycle.expRate                       = cycle.expRate .* repmat(norm,[1 size(cycle.expRate,2)]); clear norm
norm                = sum(cycle.nSpikes,2) ./ sum(sum(cycle.expRate_binned,2),3); norm(isnan(norm)) = 0;
cycle.expRate_binned                = cycle.expRate_binned .* repmat(norm,[1 size(cycle.expRate_binned,2) size(cycle.expRate_binned,3)]); clear norm
norm                = sum(cycle.nSpikes,2) ./ sum(sum(cycle.expPhase_binned,2),3); norm(isnan(norm)) = 0;
cycle.expPhase_binned               = cycle.expPhase_binned .* repmat(norm,[1 size(cycle.expPhase_binned,2) size(cycle.expPhase_binned,3)]); clear norm
norm                = sum(cycle.nSpikes,2) ./ sum(sum(cycle.expComb_binned,2),3); norm(isnan(norm)) = 0;
cycle.expComb_binned                = cycle.expComb_binned .* repmat(norm,[1 size(cycle.expComb_binned,2) size(cycle.expComb_binned,3)]); clear norm
if in.fieldRateVar>0
    norm            = sum(cycle.nSpikes,2) ./ sum(sum(cycle.expPure_binned,2),3); norm(isnan(norm)) = 0;
    cycle.expPure_binned            = cycle.expPure_binned .* repmat(norm,[1 size(cycle.expPure_binned,2) size(cycle.expPure_binned,3)]); clear norm   
    norm            = sum(cycle.nSpikes,2) ./ sum(sum(cycle.expPurC_binned,2),3); norm(isnan(norm)) = 0;
    cycle.expPurC_binned            = cycle.expPurC_binned .* repmat(norm,[1 size(cycle.expPurC_binned,2) size(cycle.expPurC_binned,3)]); clear norm
end

end



%% Function to decode movement trajectory from grid cell activity within
%  each oscillatory cycle
function[decoded]   = decodeCycles(in,pos,gridCells,cycle)

%  Extract / specify some parameters for the analysis
nCells              = in.nGCs * length(in.gridScales);
nCycles             = size(cycle.nSpikes,2);
nPhaseBins          = size(cycle.nSpikes_binned,2);
binSize             = 2;

%  Compute the mean rate function for each location bin
switch in.environment
    case {'1D','1Dlong'}
        posBins   	= 0:binSize:ceil(max(pos.x_log)/binSize)*binSize;
        [~,posBin]  = histc(pos.x_log,posBins);
        [x,y]       = ndgrid(posBin,1:nCells); clear posBin
        
        rate        = gridCells.r_rate';
        avgRate     = accumarray([x(:) y(:)],rate(:),[length(posBins)-1 nCells],@mean)'; clear rate 
        norm        = sum(cycle.nSpikes,2) ./ sum(avgRate,2);
        avgRate     = avgRate .* repmat(norm,1,size(avgRate,2)); clear norm        
        
        if in.fieldRateVar > 0            
            pure 	= gridCells.r_pure';            
            avgPure = accumarray([x(:) y(:)],pure(:),[length(posBins)-1 nCells],@mean)'; clear pure
            norm  	= sum(cycle.nSpikes,2) ./ sum(avgPure,2);
            avgPure = avgPure .* repmat(norm,1,size(avgPure,2)); clear norm
        end                
        clear x y 
    case '2D'
        maxPos      = max([pos.x_log pos.y_log],[],1);
        minPos      = min([pos.x_log pos.y_log],[],1);
        posBins     = floor(min(minPos)/binSize)*binSize : binSize : ceil(max(maxPos)/binSize)*binSize;
        posBinInd   = [pos.x_log-minPos(1) pos.y_log-minPos(2)] ./ binSize; clear minPos maxPos       
        posBinInd(posBinInd==0) = eps;
        posBinInd   = ceil(posBinInd);
        avgRate     = nan(range(posBinInd(:,1))+1,range(posBinInd(:,2))+1,nCells);        
        if in.fieldRateVar > 0
            avgPure	= nan(range(posBinInd(:,1))+1,range(posBinInd(:,2))+1,nCells);
        end
        for c       = 1 : nCells
            rate    = gridCells.r_rate(c,:)';
            avgRate(:,:,c)      = accumarray(posBinInd,rate,[],@mean); clear rate
            if in.fieldRateVar > 0
                pure 	= gridCells.r_pure(c,:)';
                avgPure(:,:,c)	= accumarray(posBinInd,pure,[],@mean); clear pure
            end
        end
        clear c posBinInd        
end

%  Decode location in each cycle / sub-cycle using firing rate across location bins
decodedLoc          = nan(nCycles,2);
binnedDecLoc        = nan(nPhaseBins,2,nCycles);
for bin             = 1 : nCycles    
    decodedLoc(bin,:)       = decodeLocation(cycle.nSpikes(:,bin),avgRate);
    if in.fieldRateVar > 0
        if bin      == 1
            decodedPure  	= nan(nCycles,2);
            binnedDecPure	= nan(nPhaseBins,2,nCycles);
        end
        decodedPure(bin,:)	= decodeLocation(cycle.nSpikes(:,bin),avgPure);
    end
    for phase       = 1 : nPhaseBins        
        binnedDecLoc(phase,:,bin)       = decodeLocation(cycle.nSpikes_binned(:,phase,bin),avgRate);
        if in.fieldRateVar > 0
            binnedDecPure(phase,:,bin)	= decodeLocation(cycle.nSpikes_binned(:,phase,bin),avgPure);
        end
    end
    clear phase
end
decoded.decLoc              = posBins(decodedLoc) + binSize/2; clear bin decodedLoc 
decoded.decLocErr           = sqrt(sum((decoded.decLoc - cycle.actLoc).^2,2)); 
decoded.binDecLoc           = posBins(binnedDecLoc) + binSize/2; clear binnedDecLoc
decoded.binDecLocErr        = squeeze(sqrt(sum((decoded.binDecLoc - cycle.actLoc_binned).^2,2)))';
decoded.binDecLocDynErr     = squeeze(decoded.binDecLoc - cycle.actLoc_binned);
if in.fieldRateVar > 0
    decoded.decPur          = posBins(decodedPure) + binSize/2; clear decodedPure
    decoded.decPurErr       = sqrt(sum((decoded.decPur - cycle.actLoc).^2,2));
    decoded.binDecPur       = posBins(binnedDecPure) + binSize/2; clear binnedDecPure
    decoded.binDecPurErr	= squeeze(sqrt(sum((decoded.binDecPur - cycle.actLoc_binned).^2,2)))';
    decoded.binDecPurDynErr	= squeeze(decoded.binDecPur - cycle.actLoc_binned);
end
clear posBins avgRate avgComb avgPure

%  Decode location in each cycle / sub-cycle using firing rate across
%  cycles, and trajectory in each cycle
cycLoc              = nan(nCycles,2);
rateLoc             = nan(nCycles,2);
phaseLoc            = nan(nCycles,2);
combLoc             = nan(nCycles,2);
decVec              = nan(nCycles,1);
if in.fieldRateVar > 0
    pureLoc         = nan(nCycles,2);
    purcLoc         = nan(nCycles,2);
end
for bin             = 1 : nCycles    
    cycLoc(bin,:)	= decodeLocation(cycle.nSpikes(:,bin),cycle.expRate);
    rate            = cycle.nSpikes_binned(:,:,bin);
    rateLoc(bin,:)	= decodeLocation(rate(:),reshape(cycle.expRate_binned,[size(cycle.expRate_binned,1)*size(cycle.expRate_binned,2) size(cycle.expRate_binned,3)]));
    phaseLoc(bin,:)	= decodeLocation(rate(:),reshape(cycle.expPhase_binned,[size(cycle.expPhase_binned,1)*size(cycle.expPhase_binned,2) size(cycle.expPhase_binned,3)]));
    combLoc(bin,:)	= decodeLocation(rate(:),reshape(cycle.expComb_binned,[size(cycle.expComb_binned,1)*size(cycle.expComb_binned,2) size(cycle.expComb_binned,3)]));
    if in.fieldRateVar > 0
        pureLoc(bin,:)	= decodeLocation(rate(:),reshape(cycle.expPure_binned,[size(cycle.expPure_binned,1)*size(cycle.expPure_binned,2) size(cycle.expPure_binned,3)]));
        purcLoc(bin,:)	= decodeLocation(rate(:),reshape(cycle.expPurC_binned,[size(cycle.expPurC_binned,1)*size(cycle.expPurC_binned,2) size(cycle.expPurC_binned,3)]));
    end
    decTraj         = decoded.binDecLoc(:,:,bin);
    bx              = regress(decTraj(:,1),[1:nPhaseBins ; ones(1,nPhaseBins)]');
    by              = regress(decTraj(:,2),[1:nPhaseBins ; ones(1,nPhaseBins)]');
    decVec(bin,1)	= atan2(by(1),bx(1)); clear bx by decTraj rate
    
end
decoded.cycLoc      = cycle.actLoc(cycLoc(:,1),:); clear bin cycLoc
decoded.cycLocErr   = sqrt(sum((decoded.cycLoc - cycle.actLoc).^2,2)); 
decoded.rateLoc     = cycle.actLoc(rateLoc(:,1),:); clear rateLoc
decoded.rateLocErr  = sqrt(sum((decoded.rateLoc - cycle.actLoc).^2,2)); 
decoded.phaseLoc    = cycle.actLoc(phaseLoc(:,1),:); clear phaseLoc 
decoded.phaseLocErr = sqrt(sum((decoded.phaseLoc - cycle.actLoc).^2,2)); 
decoded.combLoc   	= cycle.actLoc(combLoc(:,1),:); clear combLoc 
decoded.combLocErr  = sqrt(sum((decoded.combLoc - cycle.actLoc).^2,2)); 
if in.fieldRateVar > 0
    decoded.pureLoc   	= cycle.actLoc(pureLoc(:,1),:); clear pureLoc
    decoded.pureLocErr  = sqrt(sum((decoded.pureLoc - cycle.actLoc).^2,2));    
    decoded.combPureLoc     = cycle.actLoc(purcLoc(:,1),:); clear purcLoc
    decoded.combPureLocErr  = sqrt(sum((decoded.combPureLoc - cycle.actLoc).^2,2));    
end
decoded.decVec      = decVec; clear decVec
decoded.vecError    = angle(exp(1i*cycle.realVec)./exp(1i*decoded.decVec));

%  Then do linear regression on data from half of the cycles, use that to
%  predict running speed in the other half
spikesPerCycle      = sum(cycle.nSpikes)';
b                   = regress(cycle.meanSpeed(1:2:end),[spikesPerCycle(1:2:end) ones(round(nCycles/2),1)]);
decoded.predSpeed   = spikesPerCycle(2:2:end).*b(1) + b(2); clear b
decoded.speedError  = cycle.meanSpeed(2:2:end) - decoded.predSpeed; clear spikesPerCycle

%  Then do the same for 'anxiety'
b                   = regress(cycle.meanAnxiety(1:2:end),[1./cycle.cycleLength(1:2:end) ones(round(nCycles/2),1)]);
decoded.predAnx     = b(1)./cycle.cycleLength(2:2:end) + b(2); clear b
decoded.anxError    = (decoded.predAnx ./ cycle.meanAnxiety(2:2:end)) - 1;

end



%% Maximum likelihood estimation function for decoding location based on 
%  firing rates of spatial cells
function[decLoc]    = decodeLocation(firingRates,rateFunction)

expExpRate          = exp(-rateFunction);
factRateMaps        = factorial(firingRates);
if size(rateFunction,3)==1
    currK       	= repmat(firingRates,[1 size(rateFunction,2)]); clear firingRates
    fact_currK      = repmat(factRateMaps,[1 size(rateFunction,2)]); clear factRateMaps
    pval_contrib    = ((rateFunction.^currK)./fact_currK) .* expExpRate; clear expExpRate currK fact_currK rateFunction
    pval            = prod(pval_contrib); clear pval_contrib
    peaks           = find(pval==max(pval(:))); clear pval
    if ~isempty(peaks)
        decLoc     	= [peaks(randi(length(peaks),1)) 1]; 
    else
        decLoc      = nan;
    end
    clear peaks
else
    currK           = repmat(permute(firingRates,[3 2 1]), [size(rateFunction,1),size(rateFunction,2), 1]); clear firingRates
    fact_currK      = repmat(permute(factRateMaps,[3 2 1]), [size(rateFunction,1),size(rateFunction,2), 1]); clear factRateMaps
    pval_contrib    = ((rateFunction.^currK)./fact_currK) .* expExpRate; clear expExpRate currK fact_currK
    pval            = prod(pval_contrib,3); clear pval_contrib
    peaks           = find(pval==max(pval(:))); clear pval
    if ~isempty(peaks)
        [decLoc(1),decLoc(2)]	= ind2sub([size(rateFunction,1),size(rateFunction,2)],peaks(randi(length(peaks),1)));
    else
        decLoc      = [nan nan];
    end
    clear rateFunction peaks
end
end



%% Function to plot movement trajectory decoding results
function[] = plotResults(cycle,cycleDecoding)

nBins   = 30;       % Number of histogram bins
vThresh = 5;        % Speed threshold for location and direction decoding

figure
subplot(2,3,1)
xAx     = linspace(0,50,nBins);
histDat = hist(cycleDecoding.decLocErr(cycle.meanSpeed>=vThresh),xAx);
histDat = histDat./sum(histDat);
bar(xAx,histDat,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
hold on
plot([median(cycleDecoding.decLocErr(cycle.meanSpeed>=vThresh)) median(cycleDecoding.decLocErr(cycle.meanSpeed>=vThresh))],[0 1],'r','LineWidth',2)
hold off
xlabel('Absolute location decoding error (cm)','FontSize',18)
ylabel('Relative Frequency','FontSize',18)
text(10,0.75,['Median error = ' num2str(median(cycleDecoding.decLocErr(cycle.meanSpeed>=vThresh)),3) 'cm'],'FontSize',18);
axis square
clear histDat
title('Location (rate only)','FontSize',24)

subplot(2,3,2)
histDat = hist(cycleDecoding.rateLocErr(cycle.meanSpeed>=vThresh),xAx);
histDat = histDat./sum(histDat);
bar(xAx,histDat,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
hold on
plot([median(cycleDecoding.rateLocErr(cycle.meanSpeed>=vThresh)) median(cycleDecoding.rateLocErr(cycle.meanSpeed>=vThresh))],[0 1],'r','LineWidth',2)
hold off
xlabel('Absolute location decoding error (cm)','FontSize',18)
ylabel('Relative Frequency','FontSize',18)
text(10,0.75,['Median error = ' num2str(median(cycleDecoding.rateLocErr(cycle.meanSpeed>=vThresh)),3) 'cm'],'FontSize',18);
title('Location (rate only, in phase bins)','FontSize',24)
axis square
clear histDat

subplot(2,3,3)
histDat = hist(cycleDecoding.combLocErr(cycle.meanSpeed>=vThresh),xAx);
histDat = histDat./sum(histDat);
bar(xAx,histDat,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
hold on
plot([median(cycleDecoding.combLocErr(cycle.meanSpeed>=vThresh)) median(cycleDecoding.combLocErr(cycle.meanSpeed>=vThresh))],[0 1],'r','LineWidth',2)
hold off
xlabel('Absolute location decoding error (cm)','FontSize',18)
ylabel('Relative Frequency','FontSize',18)
text(10,0.75,['Median error = ' num2str(median(cycleDecoding.combLocErr(cycle.meanSpeed>=vThresh)),3) 'cm'],'FontSize',18);
title('Location (rate and phase)','FontSize',24)
axis square
clear histDat xAx

subplot(2,3,4)
xAx     = linspace(-10,10,nBins);
histDat = hist(cycleDecoding.speedError,xAx);
histDat = histDat./sum(histDat);
bar(xAx,histDat,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
xlabel('Running speed decoding error (cm)','FontSize',18)
ylabel('Relative Frequency','FontSize',18)
axis square
clear histDat xAx
ylim([0 0.2])
title('Running Speed','FontSize',24)

subplot(2,3,5)
polarhistogram(cycleDecoding.vecError(cycle.meanSpeed>=vThresh),nBins,'Normalization','Probability','FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
title('Movement Direction','FontSize',24)

subplot(2,3,6)
xAx     = linspace(-10,10,nBins);
histDat = hist(cycleDecoding.anxError*100,xAx);
histDat = histDat./sum(histDat);
bar(xAx,histDat,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
xlabel('Decoding error (%)','FontSize',18)
ylabel('Relative Frequency','FontSize',18)
axis square
ylim([0 0.15])
clear histDat xAx
title('Arbitrary Fourth Variable','FontSize',24)

end