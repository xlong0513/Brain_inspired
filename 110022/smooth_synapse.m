function ys = smooth_synapse(y,trise,tdecay)

% smooth_synapse: smoothing of spike signal with a synaptic kernel
% Convolves the vector y with a synaptic double exponential kernel 
% (1-exp(-t/trise)).*exp(-t/tdecay). This smoothing procedure keeps  the 
% causality in the reponses.
% 
% created: Jakob Heinzle 01/07

y = makecolumn(y);

% Compute synaptic kernel
 Nxs = ceil(3*tdecay);
 
 xs = [-Nxs:Nxs]';
 g  = (1-exp(-xs/trise)).*exp(-xs/tdecay);
 g(find(xs<=0))=0;
 g  = g/sum(g);            
 % Add a mirror image at beginning and end to avoid edge effects.
 lo_edge = y(Nxs:-1:1);
 hi_edge = y(end:-1:end-Nxs+1);
 y = [lo_edge; y; hi_edge];
 
 % compute the convolution of the vectors
     yst = conv(y,g);        % convolve the 'spikes' with the kernel
     Nc1 = length(yst);           
     ys = yst(2*Nxs+1:Nc1-2*Nxs+1);    % resize smoothed vector    





