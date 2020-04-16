% SETUPQUAD.M
%
% Infer full compartmental model given only access to the voltage in the
% compartments. This code is released in conjunction with the paper 
%
%	Huys QJM, Ahrens M and Paninski L (2006): Efficient estimation of
%	detailed single-neurone models
%
% and can be downloaded from 
%
%	http://www.gatsby.ucl.ac.uk/~qhuys/code.html
%
% This script uses the voltage and current traces from GETDATA.M and sets up the
% matrices and vectors for the quadratic programme to be solved in MAIN.M: 
%
%	ahat = argmin_a (1/2*a'Aa + fa) s.t. a>0
%
% This script is called by MAIN.M
%
% Copyright Quentin Huys 2006
%=========================================================================
% 			set up ahat = argmin_a (aAa + ba) s.t. a>0
%=========================================================================

%..............open probabilities given voltage...........................
fprintf('................ getting channel current shapes Jc \n')
% this just involves pasting the right bits from Jcgen into the mixing 
% matrix Jc
nch = length(g(:));
Jc = sparse(nc*T,nch);
for k=1:nc
	ind  = (k-1)*T+1:k*T;
	ind3 = (k-1)*nch/nc+1:k*nch/nc;
	ind2 = k + ([1:nch/nc]-1)*nc;
	Jc(ind,ind3) = Jcgen(:,ind2);	
end

%..............input current  ..........................................
Ji = reshape(I,nc*T,1);

%..............intercompartmental current shapes...........................
fprintf('................ getting intercompartmental current shapes \n')
Jf = sparse(nc*T,nc-1);
[chil,p] = find(tril(Wc)==1);		% get child and parent identity
for k=1:nc-1
	Wtrue(k) = W(chil(k),p(k));
	ind = (chil(k)-1)*T+1:chil(k)*T;
	Jf(ind,k) = V(:,p(k)) - V(:,chil(k))  ;
	ind = (p(k)-1)*T+1:p(k)*T;
	Jf(ind,k) = V(:,chil(k)) - V(:,p(k)) ;
end


%..............generate Hessian A and vector f....................
fprintf('................ generating Hessian and vector f\n')

J = sparse([Jc Jf Ji]);

A = 2*J'*J *delta/sigma^2;
f = - 2*J'*dv/sigma^2 ;

atrue 	= [reshape(g',nch,1); Wtrue' ; R]; % vector of true parameters
