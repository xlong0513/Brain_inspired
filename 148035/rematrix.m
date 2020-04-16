%% AUTHOR
% Praveen K. Pilly (praveen.pilly@gmail.com)
%
%% LICENSE POLICY
% Written by Praveen K. Pilly, Center for Adaptive Systems, Center for Computational Neuroscience and Neural Technology, Boston University
% Copyright 2013, Trustees of Boston University
%
% Permission to use, copy, modify, distribute, and sell this software and its documentation for any purpose is hereby granted
% without fee, provided that the above copyright notice and this permission notice appear in all copies, derivative works and
% associated documentation, and that neither the name of Boston University nor that of the author(s) be used in advertising or
% publicity pertaining to the distribution or sale of the software without specific, prior written permission. Neither Boston
% University nor its agents make any representations about the suitability of this software for any purpose. It is provided "as
% is" without warranty of any kind, either express or implied. Neither Boston University nor the author indemnify any
% infringement of copyright, patent, trademark, or trade secret resulting from the use, modification, distribution or sale of
% this software.
%
%% LAST MODIFIED
% April 5, 2013

function y=rematrix(x,I,J)
% [~] combination of repmat and permute

y=zeros([I J size(x)]);

for i=1:size(x,1)
    for j=1:size(x,2)
        y(:,:,i,j)=x(i,j);
    end
end

return