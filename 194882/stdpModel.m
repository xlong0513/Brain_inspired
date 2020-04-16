function w = stdpModel(~,w,opt)
% stdpModel
%   t       - Time in msec.
%   w       - This weight models the synaptic strength.
%   opt     - Structure with fields:
%             * TimePre  - Spiking times from the pre-synaptic cell.
%             * TimePost - Spiking times from the post-synaptic cell.
%
% RETURNS
%   w       - Adapted weight.
%
% DESCRIPTION
%   Spike timing dependent plasticity (STDP) model with synaptic 
%   potentiation and synaptic deperession.

%   Florian Raudies, 09/07/2014, Boston University.

TAU_PLUS    = 10;   % msec
TAU_MINUS   = 10;   % msec
TAU_W       = 10;   % msec
A_PLUS      = 1.2;  % Amplitude for long term potentation (LTP).
A_MINUS     = -.4;  % Amplitude for long term depression (LTD).
W_MIN       = 0;    % Minimum weight value.
W_MAX       = 1;    % Maximum weight value.
eta         = 1/(TAU_W/opt.dt); % Learning rate.

% Get time stamps of spikes of pre-synaptic and post-synaptic cells.
TimePre     = opt.TimePre;
TimePost    = opt.TimePost;
nPre        = length(TimePre);
nPost       = length(TimePost);

% Compute the time difference.
Delta       = repmat(TimePost(:),[1 nPre]) - repmat(TimePre(:)',[nPost 1]);
Pos         = Delta>0; % Pre before post.
Neg         = Delta<0; % Post before pre.

% Note that W_MAX can be overshot when having a large input signal!
w           = w + eta * ( (W_MAX-w)*A_PLUS*sum(exp(-Delta(Pos)/TAU_PLUS)) ...
                        - (W_MIN-w)*A_MINUS*sum(exp(+Delta(Neg)/TAU_MINUS)));
