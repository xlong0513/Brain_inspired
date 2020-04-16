%% AUTHOR
% Praveen K. Pilly (praveen.pilly@gmail.com)
%
%% REFERENCE
% Pilly, P. K., & Grossberg, S. (2013). Spiking neurons in a hierarchical
% self-organizing map model can learn to develop spatial and temporal
% properties of entorhinal grid cells and hippocampal place cells. PLoS One, 8(4), e60599.
%
%% LICENSE POLICY
% Written by Praveen K. Pilly, Center for Adaptive Systems, Center for Computational Neuroscience and Neural Technology, Boston University
% Copyright 2013, Trustees of Boston University
%
% Permission to use, copy, modify, distribute, and sell this software and its documentation for any purpose is hereby granted
% without fee, provided that the above copyright notice and this permission notice appear in all copies, derivative works and
% associated documentation, and that neither the name of Boston University nor that of the author(s) be used in advertising or
% publicity pertaining to the distribution or sale of the software without specific, prior written permission. Neither Boston
% University nor its agents make any representations about the suitability of this software for any purpose. It is provided "as
% is" without warranty of any kind, either express or implied. Neither Boston University nor the author indemnify any
% infringement of copyright, patent, trademark, or trade secret resulting from the use, modification, distribution or sale of
% this software.
%
%% LAST MODIFIED
% April 5, 2013

clear all
pack

datafile='sim1';
dummy_datafile=[datafile '_dummy'];

rand('state',sum(100*clock))

dt=0.002; % time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increasing the temporal resolution of the trajectory to 2 ms
% load 11207-21060501+02_t6c1.mat % TRAJECTORY 1 (~20 min); Sargolini et al. (2006)
load 11084-03020501_t2c1 % TRAJECTORY 2 (~10 min); Sargolini et al. (2006)

T=length(x1);
N_interp=10; % no change in time (20 ms to 2 ms (dt))

for j=1:T-7
    t=j+6; % discarding the first 6 location points

    Pos_x(N_interp*(j-1)+1:N_interp*(j-1)+N_interp)=x1(t)+(x1(t+1)-x1(t))*[0:(N_interp-1)]/N_interp;
    Pos_y(N_interp*(j-1)+1:N_interp*(j-1)+N_interp)=y1(t)+(y1(t+1)-y1(t))*[0:(N_interp-1)]/N_interp;

end

%--------------------------
%% NOT a necessary step if the same trajectory is presented in each trial:
%% Ensuring the first location point is the origin (0,0)
qx=Pos_x(1); qy=Pos_y(1);
v_mean=15; % average speed
kx=v_mean*qx*dt/sqrt(qx^2+qy^2); % a straight path covered at a constant speed
ky=v_mean*qy*dt/sqrt(qx^2+qy^2);
Pos_x=[kx*[0:ceil(sqrt(qx^2+qy^2)/(v_mean*dt))] Pos_x];
Pos_y=[ky*[0:ceil(sqrt(qx^2+qy^2)/(v_mean*dt))] Pos_y];

% getting rid of NaN's
Pos_x=Pos_x(~isnan(Pos_x));
Pos_y=Pos_y(~isnan(Pos_y));

clear qx qy kx ky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIMULATION SETTINGS
T=length(Pos_x);
Sf=[20 35 50]; % Input stripe spacings (in cm)
na=5; % Number of spatial phases for stripe cells
sig_g=0.07*Sf; % standard deviation of stripe subfield, based on activity drop to ~20%
StOri=-90:10:80; % Stripe directions (in deg)
Nori=length(StOri);
ns=length(Sf);
ng=3; % Number of local populations in medial entorhinal cortex (MEC)
Kg=100; % Number of map cells in each MEC population
Kp=100; % Number of map cells in the hippocampal cortex (HC) population

Nlap=1; % Number of learning trials
DiffTraj=1; % 1 if the animat traverses on a novel trajectory each trial; 0 otherwise
rot=2*pi*rand(1,Nlap); % set of angles by which the /data/ trajectory is rotated around the midpoint(/origin) to obtain a set of novel trajectories for the various trials

Wg=0.1*rand(Nori,na,ns,Kg); % Initialization of weights from stripe cells to MEC map cells
Wp=0.1*(Nori*na/(Kg*ng))*rand(Kg,ng,Kp); % Initialization of weights from MEC map cells to HC map cells
%% Scaling of initial weights from MEC to HC, so that the average net bottom-up initial weight is conserved between the two hierarchical levels
Wg_init=Wg;
Wp_init=Wp;

Box_len=100; % length of box/enclosure (in cm)
scl=2.5; % --> 2.5 X 2.5 cm spatial bins

Tmap_spc=zeros(ceil(Box_len/scl),ceil(Box_len/scl),Nlap); % Spatial occupancy map
Tmap_dir=zeros(720,Nlap); % Directional occupancy map; 360 deg divided into 720 bins

Ttr=zeros(1,Nlap); % to track time taken to run each trial

sFRmax=1000./Sf; % Scale-dependent peak firing rates of stripe cells
% Make copies for stripe directions and spatial phases
str_fr=repmat(sFRmax',[1 Nori na]);
% Rearrange dimensions to match the variable (g) for stripe cell
% activations
str_fr=permute(str_fr,[2 3 1]);

% CELLULAR PARAMETERS (for MEC and HC map cells)
Cm=0.001; % membrane capacitance (in mF/cm^2)
g_leak=0.0005; % leak channel conductance (in mS/cm^2)
g_nmda=0.025; % maximal NMDA channel conductance (in mS/cm^2)
g_gaba=0.0125; % maximal GABA channel conductance (in mS/cm^2)
E_leak=-65; % reversal potential for leak channel (in mV)
E_nmda=0; % reversal potential for NMDA channel (in mV)
E_gaba=-70; % reversal potential for GABA channel (in mV)
V_th=-50; % voltage threshold for firing a spike (in mV)
V_reset=-60; % reset potential after a spike (in mV); no absolute refractory period
V_rest=E_leak; % resting potential (in s)

tau_ampa=0.005; % time constant of AMPA channel (in s)
tau_gaba=0.01; % time constant of GABA channel (in s)
tau_nmda=0.05; % time constant of NMDA channel (in s)
tau_learning=0.05; % time constant of learning gate (in s)

lamda_w=0.001; % learning rate

% MODEL
for lap=1:Nlap

    if lap==1
        Rmaps=zeros(ceil(Box_len/scl),ceil(Box_len/scl),Nori,na,ns); % Stripe cell spatial rate maps
        Rmapg=zeros(ceil(Box_len/scl),ceil(Box_len/scl),Kg,ng); % MEC map cell spatial rate maps
        Rmapp=zeros(ceil(Box_len/scl),ceil(Box_len/scl),Kp); % HC map cell spatial rate maps
        Dmaps=zeros(720,Nori,na,ns); % Stripe cell directional rate maps
        Dmapg=zeros(720,Kg,ng); % MEC map cell directional rate maps
        Dmapp=zeros(720,Kp); % HC map cell directional rate maps
        save(dummy_datafile,'Rmaps','Rmapg','Rmapp','Dmaps','Dmapg','Dmapp')
        clear Rmaps Rmapg Rmapp Dmaps Dmapg Dmapp
    end

    Rmap_s=zeros(ceil(Box_len/scl),ceil(Box_len/scl),Nori,na,ns); % Stripe cell spatial rate maps
    Rmap_g=zeros(ceil(Box_len/scl),ceil(Box_len/scl),Kg,ng); % MEC map cell spatial rate maps
    Rmap_p=zeros(ceil(Box_len/scl),ceil(Box_len/scl),Kp); % HC map cell spatial rate maps
    Dmap_s=zeros(720,Nori,na,ns); % Stripe cell directional rate maps
    Dmap_g=zeros(720,Kg,ng); % MEC map cell directional rate maps
    Dmap_p=zeros(720,Kp); % HC map cell directional rate maps

    %---------------------
    % Spiking event registers
    Sspk=cell(Nori,na,ns); % for stripe cells
    Gspk=cell(Kg,ng); % for MEC map cells
    Pspk=cell(Kp,1); % for HC map cells
    %--------------------

    Sg1_old=zeros(Nori,na,ns); % Stripe cell-triggered AMPA channel gate
    Sg1_oldN=zeros(Nori,na,ns); % Stripe cell-triggered NMDA channel conductance
    Sg_oldL=zeros(Nori,na,ns); % Stripe cell-triggered trace for learning
    S_spk=zeros(Nori,na,ns); % = 1 or 0, depending on spiking of stripe cells

    Vg_old=V_rest*ones(Kg,ng); % Initialization of potential to resting value
    Gg1_old=zeros(Kg,ng); % MEC map cell-triggered AMPA channel gate
    Gg1_oldN=zeros(Kg,ng); % MEC map cell-triggered NMDA channel conductance
    Gg2_old=zeros(Kg,ng); % MEC map cell-triggered GABA channel conductance
    Gg_oldL=zeros(Kg,ng); % MEC map cell-triggered learning gate (also, trace for learning)
    G_spk=zeros(Kg,ng); % = 1 or 0, depending on spiking of MEC map cells

    Vp_old=V_rest*ones(Kp,1); % Initialization of potential to resting value
    Pg2_old=zeros(Kp,1); % HC map cell-triggered GABA channel conductance
    Pg_oldL=zeros(Kp,1); % HC map cell-triggered learning gate (also, trace for learning)
    P_spk=zeros(Kp,1); % = 1 or 0, depending on spiking of HC map cells

    % record of weights at the beginning of the trial
    Wg_0=Wg;
    Wp_0=Wp;

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if lap==1 | DiffTraj==1
        % Rotation of the /data/ trajectory around the midpoint(/origin)
        Posx=Pos_x*cos(rot(lap))-Pos_y*sin(rot(lap));
        Posy=Pos_x*sin(rot(lap))+Pos_y*cos(rot(lap));

        % Ensuring the trajectory doesn't go out of the 100 cm x 100 cm
        % enclosure
        Posx=max(Posx,-50+0.0001);
        Posx=min(Posx,50-0.0001);
        Posy=max(Posy,-50+0.0001);
        Posy=min(Posy,50-0.0001);

        %~~~~~~~~~~~~~~~~
        vel=sqrt( ( Posx(2:T)-Posx(1:(T-1)) ).^2+( Posy(2:T)-Posy(1:(T-1)) ).^2 )/dt; % linear speed
        thet=(atan2(Posy(2:T)-Posy(1:(T-1)),Posx(2:T)-Posx(1:(T-1)))+pi)*180/pi; % body direction
        % NOTE: Head direction is assumed to be tangential to the
        % trajectory at all times
        %~~~~~~~~~~~~~~~~

        tic
        % NOTE: Thanks to Eric Zilli (zilli@bu.edu) for the optimized code
        % below towards pre-calculating stripe cell activations (g)
        % <<< <<< <<<< <<<<<
        g=zeros(Nori,na,ns,T);
        % "Project the trajectory locations onto the preferred direction
        % vectors"
        H = [cos(deg2rad(StOri')) sin(deg2rad(StOri'))];
        projTraj = H*[Posx; Posy];
        % "Make a copy of the directional displacement for each spatial
        % phase"
        projTraj = repmat(projTraj, [1 1 na]);
        % "Rearrange dimensions to match variable g"
        projTraj = permute(projTraj, [1 3 2]);
        % <<< <<< <<<< <<<<<

        for gs=1:ns
            mus = repmat(Sf(gs)*((1:na)-1)/na, [Nori 1 T]);
            g(:,:,gs,:) = stripe_comp(mus, sig_g(gs), Sf(gs), projTraj, 1);
        end
        toc
        clear projTraj mus % to clear up some memory
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end

    tic
    for t=2:T

        if mod(t,30000)==0
            fprintf('Things are moving forward! \n')
        end

        g_old=g(:,:,:,t-1); % Input stripe cell activations at time (t-1)
        S_spk=(dt*str_fr.*g_old > rand(Nori,na,ns)); % firing rates to spikes (Poisson process)

        px=Posx(t-1);
        py=Posy(t-1);
        hd=thet(t-1);

        % Stripe cell-triggered AMPA channel gate
        Sg1_old=(1-Sg1_old).*S_spk+Sg1_old;
        Sg1=Sg1_old+dt*[-Sg1_old/tau_ampa];
        % Stripe cell-triggered NMDA channel conductance
        Sg1N=Sg1_oldN+dt*[-Sg1_oldN/tau_nmda+1000*(1-Sg1_oldN).*Sg1_old];

        % Stripe cell-triggered trace for learning
        Sg_oldL=(1-Sg_oldL).*S_spk+Sg_oldL;
        SgL=Sg_oldL+dt*[-Sg_oldL/tau_learning];

        % MEC map cell-triggered AMPA channel gate
        Gg1_old=(1-Gg1_old).*G_spk+Gg1_old;
        Gg1=Gg1_old+dt*[-Gg1_old/tau_ampa];
        % MEC map cell-triggered NMDA channel conductance
        Gg1N=Gg1_oldN+dt*[-Gg1_oldN/tau_nmda+1000*(1-Gg1_oldN).*Gg1_old];

        % MEC map cell-triggered GABA channel conductance
        Gg2_old=(1-Gg2_old).*G_spk+Gg2_old;
        Gg2=Gg2_old+dt*[-Gg2_old/tau_gaba];

        % MEC map cell-triggered learning gate
        Gg_oldL=(1-Gg_oldL).*G_spk+Gg_oldL;
        GgL=Gg_oldL+dt*[-Gg_oldL/tau_learning];

        %% Summation of bottom-up excitatory inputs to MEC map cells
        %% (limited to scale)
        S_bu=squeeze(sum( sum( Wg.*repmat(Sg1_oldN,[1 1 1 Kg]),2),1));

        %% Summation of recurrent inhibitory inputs to MEC map cells
        %% (limited to population)
        G_inh=repmat(sum(Gg2_old,1),Kg,1)-Gg2_old;

        % HC map cell-triggered GABA channel conductance
        Pg2_old=(1-Pg2_old).*P_spk+Pg2_old;
        Pg2=Pg2_old+dt*[-Pg2_old/tau_gaba];

        % HC map cell-triggered learning gate
        Pg_oldL=(1-Pg_oldL).*P_spk+Pg_oldL;
        PgL=Pg_oldL+dt*[-Pg_oldL/tau_learning];

        %% Summation of bottom-up excitatory inputs to HC map cells
        G_bu=squeeze(sum(sum(Wp.*repmat(Gg1_oldN,[1 1 Kp]),2),1));
        %% Summation of recurrent inhibitory inputs to HC map cells
        P_inh=repmat(sum(Pg2_old),1,Kp)'-Pg2_old;

        %%---------------------
        % Recording of spiking events (for the last trial)
        if lap==Nlap
            for i1=1:Nori
                for i2=1:na
                    for i3=1:ng
                        if S_spk(i1,i2,i3)==1
                            Sspk{i1,i2,i3}=[Sspk{i1,i2,i3} t];
                        end
                    end
                end
            end

            for i1=1:Kg
                for i2=1:ng
                    if G_spk(i1,i2)==1
                        Gspk{i1,i2}=[Gspk{i1,i2} t];
                    end
                end
            end

            for i1=1:Kp
                if P_spk(i1)==1
                    Pspk{i1}=[Pspk{i1} t];
                end
            end
        end
        %%---------------------

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Accumulating spikes in various spatial and directional bins
        Rmap_s(ceil((px+Box_len/2)/scl),ceil((py+Box_len/2)/scl),:,:,:)=Rmap_s(ceil((px+Box_len/2)/scl),ceil((py+Box_len/2)/scl),:,:,:)+reshape(S_spk,[1 1 Nori na ns]);
        Rmap_g(ceil((px+Box_len/2)/scl),ceil((py+Box_len/2)/scl),:,:)=Rmap_g(ceil((px+Box_len/2)/scl),ceil((py+Box_len/2)/scl),:,:)+reshape(G_spk,[1 1 Kg ng]);
        Rmap_p(ceil((px+Box_len/2)/scl),ceil((py+Box_len/2)/scl),:)=Rmap_p(ceil((px+Box_len/2)/scl),ceil((py+Box_len/2)/scl),:)+reshape(P_spk,[1 1 Kp]);

        Dmap_s(ceil(hd/0.5),:,:,:)=Dmap_s(ceil(hd/0.5),:,:,:)+reshape(S_spk,[1 Nori na ns]);
        Dmap_g(ceil(hd/0.5),:,:)=Dmap_g(ceil(hd/0.5),:,:)+reshape(G_spk,[1 Kg ng]);
        Dmap_p(ceil(hd/0.5),:)=Dmap_p(ceil(hd/0.5),:)+reshape(P_spk,[1 Kp]);

        % 2.1220e-314 * (i) is added sometimes?!
        if isreal(Rmap_s)==0 | isreal(Rmap_g)==0 | isreal(Rmap_p)==0 | isreal(Dmap_s)==0 | isreal(Dmap_g)==0 | isreal(Dmap_p)==0
            disp('Something fishy going on!')
        end

        % Tracking time spent by the animat in various bins
        Tmap_spc(ceil((px+Box_len/2)/scl),ceil((py+Box_len/2)/scl),lap)=Tmap_spc(ceil((px+Box_len/2)/scl),ceil((py+Box_len/2)/scl),lap)+dt;
        Tmap_dir(ceil(hd/0.5),lap)=Tmap_dir(ceil(hd/0.5),lap)+dt;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % MEC map cells (--> grid cells)
        Vg=Vg_old+(dt/Cm)*[g_leak*(E_leak-Vg_old)+g_nmda*S_bu'.*(3.708./(1+exp(-0.0174*Vg_old))).*(E_nmda-Vg_old)+g_gaba*G_inh.*(E_gaba-Vg_old)];
        G_spk=(Vg>=V_th);
        Vg=(V_reset-Vg).*G_spk+Vg;

        % HC map cells (--> place cells)
        Vp=Vp_old+(dt/Cm)*[g_leak*(E_leak-Vp_old)+g_nmda*G_bu.*(3.708./(1+exp(-0.0174*Vp_old))).*(E_nmda-Vp_old)+g_gaba*P_inh.*(E_gaba-Vp_old)];
        P_spk=(Vp>=V_th);
        Vp=(V_reset-Vp).*P_spk+Vp;

        Woldg=Wg;
        Woldp=Wp;

        % Learning: synaptic weight changes
        Wg=Woldg+dt*lamda_w*rematrix(Gg_oldL',Nori,na).*[repmat(Sg_oldL,[1 1 1 Kg])-Woldg.*repmat(rematrix(squeeze(sum(sum(Sg_oldL,2),1)),Nori,na),[1 1 1 Kg])];
        Wp=Woldp+dt*lamda_w*rematrix(Pg_oldL,Kg,ng).*[repmat(Gg_oldL,[1 1 Kp])-Woldp.*repmat(rematrix(squeeze(sum(sum(Gg_oldL,2),1)),Kg,ng),[1 1 Kp])];

        %% For next time step
        Sg1_old=Sg1;
        Sg1_oldN=Sg1N;
        Sg_oldL=SgL;

        Vg_old=Vg;
        Gg1_old=Gg1;
        Gg1_oldN=Gg1N;
        Gg2_old=Gg2;
        Gg_oldL=GgL;

        Vp_old=Vp;
        Pg2_old=Pg2;
        Pg_oldL=PgL;

    end
    clear g % to clear up some memory
    toc

    Ttr(lap)=toc;

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tic
    % to handle 'OUT OF MEMORY' problem
    load(dummy_datafile)
    Rmaps(:,:,:,:,:,lap)=Rmap_s;
    Rmapg(:,:,:,:,lap)=Rmap_g;
    Rmapp(:,:,:,lap)=Rmap_p;
    Dmaps(:,:,:,:,lap)=Dmap_s;
    Dmapg(:,:,:,lap)=Dmap_g;
    Dmapp(:,:,lap)=Dmap_p;
    save(dummy_datafile,'Rmaps','Rmapg','Rmapp','Dmaps','Dmapg','Dmapp')

    if lap~=Nlap
        clear Rmaps Rmapg Dmaps Dmapg Rmapp Dmapp
    end
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    toc

    save(datafile)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHED
Gker=fspecial('gaussian',[5 5],1); % Smoothing kernel

Tmap_fil=imfilter(Tmap_spc,Gker);
for k=1:Nlap
    Tmaps(:,:,:,:,:,k)=repmat(Tmap_fil(:,:,k),[1 1 Nori na ns]);
    Tmapg(:,:,:,:,k)=repmat(Tmap_fil(:,:,k),[1 1 Kg ng]);
    Tmapp(:,:,:,k)=repmat(Tmap_fil(:,:,k),[1 1 Kp]);
end

RateMaps=Rmaps;
for k=1:Nori
    k
    RateMaps(:,:,k,:,:,:)=imfilter(Rmaps(:,:,k,:,:,:),Gker)./(Tmaps(:,:,k,:,:,:)+eps);
end

RateMapg=Rmapg;
for k=1:Kg
    k
    RateMapg(:,:,k,:,:)=imfilter(Rmapg(:,:,k,:,:),Gker)./(Tmapg(:,:,k,:,:)+eps);
end

RateMapp=Rmapp;
for k=1:Kp
    k
    RateMapp(:,:,k,:)=imfilter(Rmapp(:,:,k,:),Gker)./(Tmapp(:,:,k,:)+eps);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNSMOOTHED
Tmap_fil=Tmap_spc;
for k=1:Nlap
    Tmaps(:,:,:,:,:,k)=repmat(Tmap_fil(:,:,k),[1 1 Nori na ns]);
    Tmapg(:,:,:,:,k)=repmat(Tmap_fil(:,:,k),[1 1 Kg ng]);
    Tmapp(:,:,:,k)=repmat(Tmap_fil(:,:,k),[1 1 Kp]);
end

URateMaps=Rmaps;
for k=1:Nori
    k
    URateMaps(:,:,k,:,:,:)=Rmaps(:,:,k,:,:,:)./(Tmaps(:,:,k,:,:,:)+eps);
end

URateMapg=Rmapg;
for k=1:Kg
    k
    URateMapg(:,:,k,:,:)=Rmapg(:,:,k,:,:)./(Tmapg(:,:,k,:,:)+eps);
end

URateMapp=Rmapp;
for k=1:Kp
    k
    URateMapp(:,:,k,:)=Rmapp(:,:,k,:)./(Tmapp(:,:,k,:)+eps);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHED (Directional)
Dker=fspecial('average',[1 29]);

DTmap_fil=squeeze(imfilter(Tmap_dir,Dker,'circular'));
for k=1:Nlap
    Tmaps_dir(:,:,:,:,k)=repmat(Tmap_dir(:,k),[1 Nori na ns]);
    Tmapg_dir(:,:,:,k)=repmat(Tmap_dir(:,k),[1 Kg ng]);
    Tmapp_dir(:,:,k)=repmat(Tmap_dir(:,k),[1 Kp]);
end

DateMaps=Dmaps;
for k=1:Nori
    k
    DateMaps(:,k,:,:,:)=imfilter(Dmaps(:,k,:,:,:),Dker,'circular')./(Tmaps_dir(:,k,:,:,:)+eps);
end

DateMapg=Dmapg;
for k=1:Kg
    k
    DateMapg(:,k,:,:)=imfilter(Dmapg(:,k,:,:),Dker,'circular')./(Tmapg_dir(:,k,:,:)+eps);
end

DateMapp=Dmapp;
for k=1:Kp
    k
    DateMapp(:,k,:)=imfilter(Dmapp(:,k,:),Dker,'circular')./(Tmapp_dir(:,k,:)+eps);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(datafile)

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Plots
% For example, spatial firing rate maps of MEC map cells belonging to scale
% 1
for k=1:Kg
    figure
    imagesc(RateMapg(:,:,k,1,Nlap)),axis image, axis off,colorbar
end