function layer = STDPtrain(layer,i,t)
% Modify synaptic strengths of neurons in layer i using spike timing
% dependent plasticity (STDP)
% Ignores inhibitory connections
delta = 10;
smax = 2;
firings = layer{i}.firings;
tf = t-delta;
if ~isempty(firings)
   fired = firings((firings(:,1)==tf),2);
else
   fired = [];
end;
for j=1:layer{i}.rows*layer{i}.columns
   if (ismember(j,fired))
      % Neuron j fired at time tf
      for k=1:length(layer)
         if (~isempty(layer{i}.S{k}) && ~isempty(layer{k}.firings))
            % Modify connections from layer k to neuron j
            firings2 = layer{k}.firings;
            % To be of interest, spike must arrive from layer k
            % within delta ms of time of j's firing
            d = layer{i}.delay{k};
            arr = firings2(:,1)+d(j,firings2(:,2))';
            % Incoming spikes before neuron j fired
            pre = firings2(arr<tf & arr>(tf-delta),2);
            % Incoming spikes after neuron j fired
            post = firings2(arr>tf & arr<(tf+delta),2);
            % Layer k's influence on neuron j
            s = layer{i}.S{k}(j,:);
            % Pick out connection strengths > 0
            cons = s>0;
            % Pick out the delays to firing
            ds = zeros(1,length(s));
            lf = length(firings2(:,2));
            x = firings2(:,2)';
            y = abs(tf-arr');
            for z=1:lf
               % Use the shortest delay, but ignore spikes before window of
               % interest
               if (ds(x(z))==0 || (y(z)>0 && y(z)<ds(x(z))))
                  ds(x(z)) = y(z);
               end
            end
            % Increase strengths for incoming spikes before j fired
            s(pre) = cons(pre).*(s(pre)+(smax-abs(s(pre)))*0.68.*((delta-ds(pre))/delta));
            % Reduce strengths for incoming spikes after j fired
            s(post) = cons(post).*(s(post)-(smax-abs(s(post)))*0.68.*((delta-ds(post))/delta));
            layer{i}.S{k}(j,:) = s;
         end
      end
   end
end