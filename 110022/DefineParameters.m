% DefineParameters defines some general parameters for the simulation
% of the FEF network.
%
% created: Jakob Heinzle 01/07

t=0;
dt    = 0.10;           % integration time step in ms
gfac  = 2;              % factor for weights calculated as poolsize*0.02              

%=====================================================================
% Cell parameters; here rest is V=0. Adapted from Salinas, Neural
% Computation
%=====================================================================

VeqE    = 74;       % 74 mV for AMPA
VeqI    = -10;      % from -20 (K) to 0 (Cl, shunting)
V_th    = 20;       % spike threshold
V_reset = 10;       % -60 from Troyer and Miller, 1997
taucorrE= 3.0;      % time constant for background E input
taucorrI= 3.0;      % time constant for background I input
tauME   = 20;       % in milliseconds; 20 from McCormick
tauMI   = 12;       % in milliseconds; 12 from McCormick
trefE   = 1.8;      % refractory period in ms
trefI   = 1.2;      % refractory period in ms

nfac=21;

% integration time steps
tstepEc = exp(-dt/taucorrE);
tstepEc1= 1 - tstepEc;
tstepEc2= sqrt(1 - tstepEc^2);
tstepIc = exp(-dt/taucorrI);
tstepIc1= 1 - tstepIc;
tstepIc2= sqrt(1 - tstepIc^2);

% number of iterations
iterations = ceil(tmax/dt);

% define a series of auxiliaries          
att_goal=11;
max_sac=0; 
sac_goal=0; 
tlast=0; 
fov=0; 
pos=0;
jskip=skip+1; jspkE=0; jspkI=0; % parameters used for checking for saccade and for graphics update.
recdt=-dt/2; ind=1; facE=1000/100;  facI=facE*4; % ausxiliary values used to store rates of the neurons.
