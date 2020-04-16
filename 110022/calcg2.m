function G=calcg2(G,W,spikes);

% calcg2 computes the product of a matrix W with the vector spikes, that contains
% only zeros and ones and adds it to the vector G. For large matrices
% sparse input "spikes", calcg2 is much faster than an ordinary
% multiplication. G=G+W*spikes;
% 
% created: Jakob Heinzle 01/07

ind=find(spikes);
G = G + sum(W(:,ind),2);