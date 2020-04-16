function plot_ws_firings(layer,i,tmax)
% Produces a raster plot of the firings in layer{i}
% Maps four quadrants of 16 by 16 array into a linear sequence of
% 256 neurons. Quadrants appear in the order NW, NE, SW, SE
firings = layer{i}.firings;


C1 = reshape(1:256,16,16);
C2 = zeros(16,16);
C2(1:64) = C1(9:16,1:8); % NW
C2(65:128) = C1(9:16,9:16); % NE
C2(129:192) = C1(1:8,1:8); % SW
C2(193:256) = C1(1:8,9:16); % SE
% Invert index
for k = 1:256
   C1(C2(k)) = k;
end
% Map to new neuron numbers
for k = 1:length(firings)
   firings(k,2) = C1(firings(k,2));
end
% Display graph
plot(firings(:,1),firings(:,2),'.');
xlabel('Time (ms)','FontSize',14);
ylabel('Neuron number','FontSize',14);
axis([1 tmax 1 256]);
grid on;
set(gca,'YTick',[1 64 128 192 256],'FontSize',12);
title(['Firings in area ',layer{i}.name],'FontSize',14);
set(gcf,'Position',[50,50,800,600]);