% This code is written by Claudia Clopath, Imperial College London %
% Please cite: Clopath et al. Journal of Neuroscience 2014
% "A Cerebellar Learning Model of Vestibulo-Ocular Reflex
%Adaptation in Wild-Type and Mutant Mice"

% This code calls VOR.m

clear all
% Parameters
N_inp = 100;                        % number of Granule cells
BH = 2.85;                          % upper bound of plasticity at the Granule Cells to Purkinje
BL = 0.85;                          % lower bound of plasticity at the Granule Cells to Purkinje
CF_noise = 3.5;                     % noise in Climbing Fibers (CF)
cf_vest = 0.03;                     % vestibular component in Climbing Fibers (CF)
w_IP = 1;                           % fixed weight from interneuron to Purkinje cells
w_MD = 0.88;                        % initial Mossy Fiber to MVN weight
win = 1.85/N_inp;                   % initial weight from Granule Cells to Purkinje
w_GP = win*ones(N_inp,1);           % weight from Granule Cells to Purkinje
delay = 100;                        % delay of the Climbing Fibers [ms]

% Simulation parameters
T_pat = floor(1000/0.6);            % Time of a cycles [ms]
day = 50;                           % time of day training [unit of cycles]
nit = 1440;                         % time of night [unit of cycles]
Ntot = day*5+nit*5;                 % total simulating time
D_P = zeros(Ntot);                  % recording phase of the eye movement
D_G = zeros(Ntot);                  % recording gain of the eye movement
previous_t = 1;                     % start simulating time at 1

% Learning rates
alphai = 3.5000e-07;                % learning rate for w_GP learning term
alphaf = 7.4697e-05;                % learning rate for w_GP forgetting term 
alphad = 5.6022e-06;                % learning rate for w_MD
    
%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate neuronal firing rates %%%%%%%%%%%
% Granule Cells G
G = zeros(N_inp, T_pat);
dist = (0.1886*cos((1+N_inp:2*N_inp)*2*pi/N_inp))+(1:N_inp)*2*pi/N_inp; % over distribution of Granule cells firing at certain phases
for i = 1:N_inp
    G(i,:)= (cos((T_pat+1:2*T_pat)*2*pi/T_pat-dist(i))+1);
end

% Interneuron In
In = 2.5/N_inp *sum(G);
In = In -mean(In)+0.85;

% Compute baseline firing rate of Purkinje cells, P
Pmean = win*ones(N_inp,1)'* G - w_IP'* In;

% Mossy Fibers, M
Mmean = 0.25;   % mean firing rate of M
M = 0.25*cos((1+T_pat*3/4:T_pat+T_pat*3/4)*2*pi/T_pat)+Mmean;

%%%%%%%%%%%%%%%%%%%%%%%%%%% START SIMULATION OVER  TIME %%%%%%%%%%%%%%%%
% Initial phase before training
% Day before train
light = 1;      % in the light condition
gain =1;        % target gain 
Simul_t = day;  % simulating time
VOR             % call the function VOR that loop over time

% Night before training
light = 0;      % in the dark
gain =1;
Simul_t = 2*nit;
VOR

% Start of the training
% Day 1
light = 1;
gain = 0;
Simul_t = day;
VOR

%Night 1
light = 0;
gain =1;
Simul_t = nit;
VOR

% Day 2
light = 1;
gain =-0.5;
Simul_t = day;
VOR

% Night 2
light = 0;
gain =1;
Simul_t = nit;
VOR

% Day 3
light = 1;
gain =-1;
Simul_t = day;
VOR

% Night 3
light = 0;
gain =1;
Simul_t = nit;
VOR

% Day 4
light = 1;
gain =-1;
Simul_t = day;
VOR

%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renormalise phase and gain
D_G(1,:)= D_G(1,:)/D_G(1,2*nit);
D_P = (D_P /T_pat)*360;
D_P = mod(360-D_P -269, 360);
D_P = mod(D_P +10, 360)-10;
% Plot gain
start_train = day+2*nit;
figure;hold on
plot((1:day),D_G(start_train+(1:day))','g',  'linewidth',2)
plot(((day+1):(2*day)), D_G(start_train+((day+nit+1):(2*day + nit)))', 'g', 'linewidth',2)
plot(((2*day+1):(3*day)), D_G(start_train+((2*day+2*nit+1):(3*day + 2*nit)))','g',  'linewidth',2)
plot(((3*day+1):(4*day)), D_G(start_train+((3*day+3*nit+1):(4*day + 3*nit)))'','g',  'linewidth',2)
plot(day+(0:0.1:1.4)*0,0:0.1:1.4, '--k')
plot(2*day+(0:0.1:1.4)*0,0:0.1:1.4, '--k')
plot(3*day+(0:0.1:1.4)*0,0:0.1:1.4, '--k')
xlim([-1 201])
ylabel('gain','fontsize',20);
xlabel('time [min]','fontsize',20);
title('day 1       day 2       day 3       day 4      ','fontsize',20);
set(gca,'fontsize',20);
% Plot phase
figure; hold on;
plot((1:day),D_P(start_train+(1:day))','g',  'linewidth',2)
plot(((day+1):(2*day)), D_P(start_train+((day+nit+1):(2*day + nit)))','g',  'linewidth',2)
plot(((2*day+1):(3*day)), D_P(start_train+((2*day+2*nit+1):(3*day + 2*nit)))', 'g', 'linewidth',2)
plot(((3*day+1):(4*day)), D_P(start_train+((3*day+3*nit+1):(4*day + 3*nit)))'','g',  'linewidth',2)
plot(day+(-5:355)*0,-5:355, '--k')
plot(2*day+(-5:355)*0,-5:355, '--k')
plot(3*day+(-5:355)*0,-5:355, '--k')
xlim([-1 201])
ylim([-5 180])
ylabel('phase [deg]','fontsize',20);
xlabel('time  [min]','fontsize',20);
title('day 1       day 2       day 3       day 4      ','fontsize',20);
set(gca,'fontsize',20);