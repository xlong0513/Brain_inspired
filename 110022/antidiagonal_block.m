function W=antidiagonal_block(n1,n2,ngroups)

% W = antidiagonal_block(n1,n2,ngroups) creates a n2xn1 matrix W 
% that consists of ngroups blocks of ones on the antidiagonal.
% 
% created: Jakob Heinzle 01/07

W=zeros(n2,n1);
bs1=n1/ngroups;
bs2=n2/ngroups;
for i=1:ngroups
    W((i-1)*bs2+1:i*bs2,(ngroups-i)*bs1+1:(ngroups-i+1)*bs1)=ones(bs2,bs1);
end
    