% update_population updates the all populations. The code is adapted from 
% (Salinas, NeuralComputation,15,1439ff). 
%
% created: Jakob Heinzle 01/07

for k=1:pops.npops
   % update input conductances first
   GE=zeros(pops.population{k}.n_neurons,1);
   GI=zeros(pops.population{k}.n_neurons,1);
   INP=zeros(pops.population{k}.n_neurons,1);
   INPNOISE=zeros(pops.population{k}.n_neurons,1);
   for uu=pops.population{k}.input_exc
      GE=GE+cons.connection{uu}.G;
   end
   for uu=pops.population{k}.input_inh
      GI=GI+cons.connection{uu}.G;
   end
   % update external inputs
   for uu=pops.population{k}.input_external
      if t>inputs.external{uu}.t_off
         %dont add external input.
      elseif t>inputs.external{uu}.t_trans_off
         INP=INP+inputs.external{uu}.sustained_level*inputs.external{uu}.ExtInp;
         INPNOISE=INPNOISE+inputs.external{uu}.sustained_level*inputs.external{uu}.NoiseExtInp.*randn(pops.population{k}.n_neurons,1);
      elseif t>inputs.external{uu}.t_on
         INP=INP+inputs.external{uu}.ExtInp;
         INPNOISE=INPNOISE+inputs.external{uu}.NoiseExtInp.*randn(pops.population{k}.n_neurons,1);
      end
   end
   % calculate the total conductances of the neurons, the membrane voltage
   % and the spikes for excitatory neurons.
   if pops.population{k}.type=='exc';
      pops.population{k}.inEAux = pops.population{k}.inEAux*tstepEc + ...
         (pops.population{k}.bgExc + INP)*tstepEc1+ ...
         (pops.population{k}.NoiseExc*randn(pops.population{k}.n_neurons,1)+ INPNOISE)*tstepEc2;
      pops.population{k}.inIAux = pops.population{k}.inIAux*tstepIc + ...
         (pops.population{k}.bgInh)*tstepIc1 + ...
         (pops.population{k}.NoiseInh*randn(pops.population{k}.n_neurons,1))*tstepIc2;
      temp = 1./(1 + GE + GI + pops.population{k}.inEAux + pops.population{k}.inIAux);
      tauV = tauME*temp;
      Vinf = ((GE + pops.population{k}.inEAux)*VeqE + (GI + pops.population{k}.inIAux)*VeqI).*temp;
      % Calculate new values (!! refractory period !!)
      pops.population{k}.refrac = rectify(pops.population{k}.refrac - dt);
      iiref = (pops.population{k}.refrac == 0);
      pops.population{k}.Vm(iiref) = Vinf(iiref) + ...
         (pops.population{k}.Vm(iiref) - Vinf(iiref)).*exp(-dt./tauV(iiref));
      % look for neurons who crossed threshold.
      pops.population{k}.spikes = pops.population{k}.Vm > V_th;
      pops.population{k}.refrac(pops.population{k}.spikes) = trefE;
      pops.population{k}.Vm(pops.population{k}.spikes) = V_reset;
   end

   % calculate the total conductances of the neurons, the membrane voltage
   % and the spikes for inhibitory neurons.
   if pops.population{k}.type=='inh';
      pops.population{k}.inEAux = pops.population{k}.inEAux*tstepEc + ...
         (pops.population{k}.bgExc + INP)*tstepEc1 + ...
         (pops.population{k}.NoiseExc*randn(pops.population{k}.n_neurons,1)+INPNOISE)*tstepEc2;
      pops.population{k}.inIAux = pops.population{k}.inIAux*tstepIc + ...
         (pops.population{k}.bgInh)*tstepIc1 + ...
         (pops.population{k}.NoiseInh*randn(pops.population{k}.n_neurons,1))*tstepIc2;
      temp = 1./(1 + GE + GI + pops.population{k}.inEAux + pops.population{k}.inIAux);
      tauV = tauMI*temp;
      Vinf = ((GE + pops.population{k}.inEAux)*VeqE + (GI + pops.population{k}.inIAux)*VeqI).*temp;
      % Calculate new values (!! refractory period !!)
      pops.population{k}.refrac = rectify(pops.population{k}.refrac - dt);
      iiref = (pops.population{k}.refrac == 0);
      pops.population{k}.Vm(iiref) = Vinf(iiref) + ...
         (pops.population{k}.Vm(iiref) - Vinf(iiref)).*exp(-dt./tauV(iiref));
      % look for neurons who crossed threshold.
      pops.population{k}.spikes = pops.population{k}.Vm > V_th;
      pops.population{k}.refrac(pops.population{k}.spikes) = trefE;
      pops.population{k}.Vm(pops.population{k}.spikes) = V_reset;
   end
end


