function count = count_firings(layer,i,tmin,tmax,a,b)
% Count neuron firings for layer i between times tmin and tmax, for
% neurons numbered from a to b

firings = layer{i}.firings;
% Convert indices
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
F = firings(:,1) >= tmin & firings(:,1) < tmax & ...
   firings(:,2) >= a & firings(:,2) < b;
count = sum(F);