% PARAM.M
%
% Infer full compartmental model given only access to the voltage in the
% compartments. This code is released in conjunction with the paper 
%
%	Huys QJM, Ahrens M and Paninski L (2006): Efficient estimation of
%	detailed single-neurone models
%
% and can be downloaded from 
%
%	http://www.gatsby.ucl.ac.uk/~qhuys/code.html
%
% This script sets all the variables for the simulations. It is called by MAIN.M
%
% Copyright Quentin Huys 2006

%==========================================================================================
%			SIMULATION PARAMETERS
%==========================================================================================
nc = 30;			% number of compartments
delta = .01;       		% timestep in [ms] best <.05
T = 2000;			% number of time points
Tmax = T*delta;			% resulting length of recording

sigma = 160;			% SD of current noise [mA/cm^2]

Imean = -10;			% input current to soma: sinusoidal of amplitude Ivar [mA/cm^2]
Ivar = 1000;			% plus constant of amplitue Imean [mA/cm^2]
Ifreq = .1;			% frequency of stimulating (sinusoid)^2 [kHz]
noinputind = [2:nc];		% which compartments should NOT get sinusoidal current

Wscale = 200;                  	% [mS/cm^2] this implies a compartment length 
				% L = sqrt(10^5/(2*Wscale)) um for radius=1um, R_ax=1kOhm*mm 
L = sqrt(10^5/(2*Wscale));      fprintf('intercompartmental conductance implies L=%gum; 1um radius\n',L)

usederiv = 1;			% set this to zero to use the true current, rather than derivative

%==========================================================================================
%			CHANNELS AND MORPHOLOGY
%==========================================================================================
g = repmat([140,36,3],nc,1);     	% Average conductances for HH channels [mS/cm^2]
g(:,1) = g(:,1) + randn(nc,1)*10;
g(:,2) = g(:,2) + randn(nc,1)*3;
g(:,3) = g(:,3) + randn(nc,1)*.5;

R = 1;					% resistance to input current
C = 1;     				% capacitance [muF/cm^2], all parameters relative to C

inivars = [	-60	*ones(1,nc); 	% initial variables for differential equation solver
		0.0785	*ones(1,nc); 
		0.5637	*ones(1,nc);
		0.4147	*ones(1,nc)];	

ek=-77;el=-54;ena=50;			% reversal potentials
E = [ena, ek, el]; 

%==========================================================================================
%			OPTIMIZER
%==========================================================================================
% Use Matlab's own quadratic minimizer 'quadprog'. 
% This needs the optimization toolbox. 

minimizer = 'qp';

% Use minq5 to solve the quadratic programme
% minq5 can be downloaded for free from http://www.mat.univie.ac.at/~neum/software/minq/
% Results tend to be rather noisy. 
%
% minimizer = 'minq';

%==========================================================================================
%			PLOTS
%==========================================================================================

doplotcell = 0; 	% set this to 1 to plot the cell. This needs GraphViz2Mat on unix system

