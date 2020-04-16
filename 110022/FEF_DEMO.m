% FEF_DEMO simulates the FEF network. All changes to the input of the network
% need to be made in the file SetupDemo, where the task is defined.
%
% created: Jakob Heinzle 01/07

% reset the random generator of Matlab.
clear all;
rand('state',sum(100*clock));

% Define task and set some parameters and the readout variables.
SetupDemo;
DefineParameters;
ReadOutVariables;

% create the network and define populations, connections and inputs.
initialize_network;
makeall_populations;
makeall_connections;
makeall_input;

%initialize a series of counters counting the total amount
%of spikes for each population.
countE4=0;countI4=0;countE23=0;countI23=0;countE5R=0;countI5R=0;countE5B=0;countI5B=0;
countE6A=0;countE6S=0;countI6=0;countIFIX=0;countE6R=0;

% run the simulation
disp('Simulation is running.');
tic;
for j=1:iterations
   t=t+dt;
   recdt=recdt+dt;
   update_conductance;
   update_population;

   % check if there is a  threshold crossing of layer 5
   if recdt>1
      recdt=dt/2;
      ind=ind+1;
      if ind>10
         % check, where the attention is (for display)
         [max_att,att_goal]=max(mean(E23_HZ(:,ind-10:ind),2));
         if max_att<30
            att_goal=-20;
         end
         ATT(ind)=att_goal;
         % check for saccades.
         [max_sac,sac_goal]=max(mean(E5B_HZ(:,ind-10:ind),2)); % take the mean firing over 10 ms.
      end
      if (max_sac>50) % threshold crossing ?
         if (t-tlast>50) % minimum interval between saccades, to not detect the same saccade several times.
            sac_shift=round(sac_goal-(nfac+1)/2);
            disp(['saccade to position',int2str(sac_shift),' at time ',int2str(round(t))]);
            fov=fov+sac_shift; % calculate the new position of the fovea.
            tlast=t;
            update_input_saccade; % update input according to saccade.
         end
      end
   end

   % save the spike pattern for populations in the FEF
   E4(countE4+1:countE4+length(find(pops.population{1}.spikes)),:)=[t*ones(length(find(pops.population{1}.spikes)),1) find(pops.population{1}.spikes)]; countE4=countE4+length(find(pops.population{1}.spikes));
   I4(countI4+1:countI4+length(find(pops.population{2}.spikes)),:)=[t*ones(length(find(pops.population{2}.spikes)),1) find(pops.population{2}.spikes)]; countI4=countI4+length(find(pops.population{2}.spikes));
   E23(countE23+1:countE23+length(find(pops.population{3}.spikes)),:)=[t*ones(length(find(pops.population{3}.spikes)),1) find(pops.population{3}.spikes)]; countE23=countE23+length(find(pops.population{3}.spikes));
   I23(countI23+1:countI23+length(find(pops.population{4}.spikes)),:)=[t*ones(length(find(pops.population{4}.spikes)),1) find(pops.population{4}.spikes)]; countI23=countI23+length(find(pops.population{4}.spikes));
   E5R(countE5R+1:countE5R+length(find(pops.population{5}.spikes)),:)=[t*ones(length(find(pops.population{5}.spikes)),1) find(pops.population{5}.spikes)]; countE5R=countE5R+length(find(pops.population{5}.spikes));
   I5R(countI5R+1:countI5R+length(find(pops.population{6}.spikes)),:)=[t*ones(length(find(pops.population{6}.spikes)),1) find(pops.population{6}.spikes)]; countI5R=countI5R+length(find(pops.population{6}.spikes));
   E5B(countE5B+1:countE5B+length(find(pops.population{7}.spikes)),:)=[t*ones(length(find(pops.population{7}.spikes)),1) find(pops.population{7}.spikes)]; countE5B=countE5B+length(find(pops.population{7}.spikes));
   I5B(countI5B+1:countI5B+length(find(pops.population{8}.spikes)),:)=[t*ones(length(find(pops.population{8}.spikes)),1) find(pops.population{8}.spikes)]; countI5B=countI5B+length(find(pops.population{8}.spikes));
   E6A(countE6A+1:countE6A+length(find(pops.population{10}.spikes)),:)=[t*ones(length(find(pops.population{10}.spikes)),1) find(pops.population{10}.spikes)]; countE6A=countE6A+length(find(pops.population{10}.spikes));
   E6S(countE6S+1:countE6S+length(find(pops.population{9}.spikes)),:)=[t*ones(length(find(pops.population{9}.spikes)),1) find(pops.population{9}.spikes)]; countE6S=countE6S+length(find(pops.population{9}.spikes));
   IFIX(countIFIX+1:countIFIX+length(find(pops.population{11}.spikes)),:)=[t*ones(length(find(pops.population{11}.spikes)),1) find(pops.population{11}.spikes)]; countIFIX=countIFIX+length(find(pops.population{11}.spikes));
   FOVEA(ind)=fov;
   % save the average firing rate for populatiions in the REC module
   EFp(:,ind)=EFp(:,ind)+sum(reshape(pops.population{13}.spikes,pops.population{13}.poolsize,nfac))'*facE;
   EFf(:,ind)=EFf(:,ind)+sum(reshape(pops.population{14}.spikes,pops.population{14}.poolsize,nfac))'*facE;
   EFa(:,ind)=EFa(:,ind)+sum(reshape(pops.population{15}.spikes,pops.population{15}.poolsize,nfac))'*facE;
   IF(:,ind)=IF(:,ind)+sum(reshape(pops.population{12}.spikes,pops.population{12}.poolsize,nfac))'*facI;
   ERr(:,ind)=ERr(:,ind)+sum(reshape(pops.population{17}.spikes,pops.population{17}.poolsize,nfac))'*facE;
   ERb(:,ind)=ERb(:,ind)+sum(reshape(pops.population{18}.spikes,pops.population{18}.poolsize,nfac))'*facE;
   IRb(:,ind)=IRb(:,ind)+sum(reshape(pops.population{19}.spikes,pops.population{19}.poolsize,nfac))'*facI;

   E4_HZ(:,ind)=E4_HZ(:,ind)+sum(reshape(pops.population{1}.spikes,pops.population{1}.poolsize,nfac))'*facE;
   E23_HZ(:,ind)=E23_HZ(:,ind)+sum(reshape(pops.population{3}.spikes,pops.population{3}.poolsize,nfac))'*facE;
   E5B_HZ(:,ind)=E5B_HZ(:,ind)+sum(reshape(pops.population{7}.spikes,pops.population{7}.poolsize,nfac))'*facE*2.5;
   E5R_HZ(:,ind)=E5R_HZ(:,ind)+sum(reshape(pops.population{5}.spikes,pops.population{5}.poolsize,nfac))'*facE*2.5;

   % draw the population rates of layers 4, 2/3 and 5 and the input.
   jskip=jskip+1;
   if jskip>skip/dt, jskip=0.5; 
      if ind>10
         if t>inputs.external{1}.t_off;
            pp=0;
         elseif t>inputs.external{1}.t_trans_off;
            pp=inputs.external{1}.sustained_level;
         elseif t>inputs.external{1}.t_on;
            pp=1;
         else
            pp=0;
         end
         set(hInp,'YData',inputs.external{1}.ExtInp(10:100:2040)/inputs.external{1}.MeanInp*pp);
         set(hL4,'YData',mean(E4_HZ(:,ind-5:ind),2));
         set(hL23,'YData',mean(E23_HZ(:,ind-5:ind),2));
         set(hL5,'YData',mean(E5B_HZ(:,ind-5:ind)+E5R_HZ(:,ind-5:ind),2));

         if att_goal>-15
            set(hText,'String',[int2str(round(t)),' ms; Attention on position',int2str(att_goal-11)]);
         else
            set(hText,'String',[int2str(round(t)),' ms; No attention allocated']);
         end
         drawnow;
      end

   end


   if tmax<t
      break;
   end
end

% cut all the zeros at the end
E4=E4(1:max(find(E4(:,1))),:);
I4=I4(1:max(find(I4(:,1))),:);
E23=E23(1:max(find(E23(:,1))),:);
I23=I23(1:max(find(I23(:,1))),:);
E5R=E5R(1:max(find(E5R(:,1))),:);
I5R=I5R(1:max(find(I5R(:,1))),:);
E5B=E5B(1:max(find(E5B(:,1))),:);
I5B=I5B(1:max(find(I5B(:,1))),:);
E6A=E6A(1:max(find(E6A(:,1))),:);
E6S=E6S(1:max(find(E6S(:,1))),:);
IFIX=IFIX(1:max(find(IFIX(:,1))),:);
simtime=toc;
eval(savepath);

% calculate population rates and print the results.
df=diff(FOVEA);
sac_positions=df(find(df));
spike_to_rate;
plot_popfov;
for i=1:length(print_positions)
   pp=print_positions(i)+11;
   plot_popE;
   plot_popI;
end
