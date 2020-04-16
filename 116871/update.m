function layer = update(layer,i,t,Dmax)
% Update membrane potential v and reset rate u for neurons in layer i
% using Izhikevich's neuron model
tau = 0.2; % Simulation time increment
% Calculate current from incoming spikes
for j=1:length(layer)
   S = layer{i}.S{j};
   if ~isempty(S)
      firings = layer{j}.firings;
      if ~isempty(firings)
         % Find incoming spikes (taking account of propagation delays)
         delay = layer{i}.delay{j};
         F = layer{i}.factor{j};
         % Sum current from incoming spikes
         k = size(firings,1);
         while (k>0 && firings(k,1)>t-Dmax)
            spikes = (delay(:,firings(k,2))==t-firings(k,1));   
            layer{i}.I(spikes) = layer{i}.I(spikes)+S(spikes,firings(k,2))*F;
            k = k-1;
         end;
         % Don't let I go below zero
         layer{i}.I = layer{i}.I.*(layer{i}.I > 0);
      end
   end
end
% Update v and u using Izhikevich's model in increments of tau
for k=1:1/tau
   v = layer{i}.v;
   u = layer{i}.u;
   layer{i}.v = v+(tau*(0.04*v.^2+5*v+140-u+layer{i}.I));
   layer{i}.u = u+(tau*(layer{i}.a.*(layer{i}.b.*layer{i}.v-u)));
   % Reset neurons that have spiked
   fired = find(layer{i}.v>=30); % indices of spikes
   if ~isempty(fired)
      layer{i}.firings = [layer{i}.firings ; t+0*fired, fired];
      layer{i}.v(fired) = layer{i}.c(fired);
      layer{i}.u(fired) = layer{i}.u(fired)+layer{i}.d(fired);
   end
end