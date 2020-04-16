function [t,r_g] = filterMyInput(lgn_struct,in_struct)

% FILTERMYINPUT - Function to evaluate the effective input to the LGN relay
% cell. The input should either be: 1. The ganglion cell firing rate, i.e.
%                                   the light stimulus filtered through 
%                                   the retinal circuit
%                                   2. The light stimulus + impulse
%                                   response of the (linear) retinal
%                                   circuit.

% unpack the parameters needed for the feedforward kernel and input
% parameters
eta_ffi = lgn_struct.eta_ffi;
tau_rg = lgn_struct.tau_rg;
tau_rig = lgn_struct.tau_rig;
Delta_rig = lgn_struct.Delta_rig;
form = in_struct.form;

% create feedforward kernel with respective time vector
tau_max = max([tau_rg tau_rig]);
n_ff = 1000;
t_ff_max = Delta_rig+4*tau_max;
dt = t_ff_max/(n_ff-1);

t_ff = 0:dt:t_ff_max;                                          % time vector
h_rg = 1/tau_rg*exp(-t_ff/tau_rg);                          % excitatory kernel
h_ffi = 1/tau_rig*exp(-(t_ff-Delta_rig)/tau_rig) .* (t_ff>=Delta_rig); % inhibitory kernel
h_ff=(h_rg - eta_ffi*h_ffi);

if strcmp(form,'vec')
    % unpack the input data and interpolate the vectors
    t_tmp = in_struct.t_in;
    r_tmp = in_struct.r_in;
    t_in=t_tmp(1):dt:t_tmp(end);
    r_in=interp1(t_tmp,r_tmp,t_in,'linear');
    clear t_tmp r_tmp
elseif strcmp(form,'fnc')
    h = in_struct.h;
    tstop = in_struct.tstop;
else
end

% 
r_g = (dt * conv(h_ff,r_in))';
t_fin=t_in(end)-t_in(1)+t_ff_max;
t=(0:dt:t_fin)';

