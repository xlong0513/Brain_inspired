function M=select_sparse_connect(M_IN,p);

% M = select_sparse_connect(M_IN,p) randomly leaves a percentage p of
% the original Matrix M_IN intact and resets all the others to zero. 
% The total strength of the connection is conserved by multiplying 
% M by 1/p,. p i s a value between 0 and 1.
%
% created: Jakob Heinzle 01/07

randhelp=rand(size(M_IN));
randhelp=randhelp<p;

M=(randhelp.*M_IN)/p;

