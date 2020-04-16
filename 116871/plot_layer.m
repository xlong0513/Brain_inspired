function plot_layer(layer,i,t)
% Displays neuron firings in layer i at time t
if ~isempty(layer{i}.firings)
   fired = layer{i}.firings(find(layer{i}.firings(:,1)==t),2);
else
   fired = [];
end;
subplot(5,5,layer{i}.pos);
cla;
s = ceil(sqrt(layer{i}.rows*layer{i}.columns));
plot(floor((fired-1)./s)+1,mod(fired-1,s)+1,'.');
title(layer{i}.name);
axis([1 s 1 s]);
axis square;
end