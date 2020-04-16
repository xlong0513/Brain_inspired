% This code is written by Claudia Clopath, Imperial College London %
% Please cite: Clopath et al. Journal of Neuroscience 2014
% "A Cerebellar Learning Model of Vestibulo-Ocular Reflex
% Adaptation in Wild-Type and Mutant Mice"

% Loop over time
Dmean = 2.25;
Dt = gain*0.25*cos((1+T_pat*3/4:T_pat+T_pat*3/4)*2*pi/T_pat)+1; % target D (ie. target MVN output)
for t = previous_t:Simul_t+ previous_t
    P = w_GP' * G - w_IP' * In;             % Purkinje cells, In is the interneuron, G the granuel cells
    D = w_MD* 2*(M-Mmean)+Dmean -M -P;      % Medial Vestibular Nuclei (MVN) cells 
    Cf = light*(Dt-D)+(M-Mmean)*cf_vest;    % Climbing Fibers (CF)
    Cf = circshift(Cf,[0,delay]);           % Delay in the CF
    
    % plasticity of G to P synapses
    w_GP = w_GP + alphai * (-1)*sum(((ones(N_inp,1)*Cf)+CF_noise*randn(N_inp, T_pat)).* G,2); % update 
    w_GP = (w_GP-(BL/N_inp)).*((w_GP-(BL/N_inp))>0) +(BL/N_inp);    % lower bound on the weights
    w_GP = (w_GP-(BH/N_inp)).*((w_GP-(BH/N_inp))<0)+(BH/N_inp);     % upper bound on the weights
    w_GP =  w_GP + alphaf *(win-w_GP);                              % update of decay
    
    % plasticity of MF to MVM synapses
    w_MD = w_MD + alphad*sum((-M+Mmean).*(P-Pmean));
    w_MD = w_MD*(w_MD>0);                   % lower bound at zero
    
    % record phase and gain of D
    [D_G( t), D_P(t)]=max(D); 
    D_G(t) = D_G(t)-mean(D);
end
previous_t = t; 