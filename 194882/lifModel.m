function V = lifModel(~, V,opt)
% lifModel
%   t   - Time in msec.
%   V   - Membrane potential in Volts.
%   opt - Structure with fields:
%         * dt      - Time simulated in one step in milli seconds.
%         * V_PEAK  - Peak value for membrane potential in Volts.
%         * V_RESET - Reset value for membrane potential in Volts.
%         * V_TH    - Threshold voltage for spike in Volts.
%
% RETURN
%   V   - Updated membrane potential in Volts.
%
% DESCRIPTION
%   Leaky Integrate and Fire (LIF) model.

%   Florian Raudies, 09/07/2014, Boston University.

dt          = opt.dt;       % Step width in milli seconds
V_PEAK      = opt.V_PEAK;   % Peak value in Volts.
V_RESET     = opt.V_RESET;  % Reset value in Volts.
V_TH        = opt.V_TH;     % Threshold value in Volts.
C_MEM       = 5.5;          % Membrane capacitance in nano Farad.
G_LEAK      = 10;           % Membrance leaky conductance in nano Siemens.

% Does any of the potentials have the peak value or above?
AbovePeak   = V>=V_PEAK;

% If they have they have the peak value or above then reset.
V(AbovePeak)= V_RESET;

% Does any of the potentials reach a value above threshold?
AboveTh     = V>V_TH;

% If a potental is above threshold set it to the peak value.
V(AboveTh)  = V_PEAK;

% Update all membrane potentials not updated thus far by using the lif eq.
Update      = ~AbovePeak & ~AboveTh;
V(Update)   = V(Update) + dt*10^-3*(G_LEAK/C_MEM*(V_RESET-V(Update)) ...
                                  + opt.I(Update)/C_MEM);
