function MR=distribute_random(M,distribution,sigma_dist);

% MR = distribute_random(M,distribution,sigma_dist) takes the nonzero values in
% M as the mean of a distribution and replaces its value with one drawn
% from the distribution. 
% distribution specifies the kind of distribution,
% that is used. 0: no distribution, 1: uniform (sigma_dist*mean is the
% range), 2: gaussian (sigma_dist*mean gives the HWHH of the distribution),
% the weights are normalized, to keep the sum of the matrix constant.
%
% created: Jakob Heinzle 01/07

if distribution==0
    MR=M;
    return;
end

if distribution==1
MR=M.*(1+2*sigma_dist*(rand(size(M))-0.5));
return;
end

if distribution==2
MR=M.*(1+normrnd(0,sigma_dist,size(M)));
MR=MR.*(MR>0);
MR=MR/sum(MR(:))*sum(M(:));
return;
end