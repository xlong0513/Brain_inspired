% GETDATA.M
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
% This script uses parameters in PARAM.M and is called by MAIN.M. It 
%
%	a) generates the dendritic tree
%	b) simulates the cell
%	c) collects current and voltage traces needed for inference
%
% Copyright Quentin Huys 2006

%-------------------------------------------------------------------------
% 			GENERATE DENDRITIC TREE
%-------------------------------------------------------------------------
fprintf('................ generating dendritic tree \n')
Wc = zeros(nc);
W = zeros(nc);
% iterate through the compartments, for each find a parent to attach it to
for k=2:nc			
	rp = rand; 
	if rp < .5; 	par(k) = ceil(rand(1)*(k-1));
	else 		par(k) = k-1;
	end
	% connectivity matrix
	Wc(k,par(k)) = 1;
	% generate random conductance
	W(k,par(k)) = Wscale + randn(1)*Wscale/10; 
end
Wc = (Wc+Wc');		% the connectivity matrix needs to be symmetric 
W = (W+W');

%-------------------------------------------------------------------------
% 			NOISE AND INPUT CURRENT
%-------------------------------------------------------------------------
Ipar{1}=Imean; Ipar{2}=Ivar; Ipar{3}=Ifreq; Ipar{4}=noinputind;
for t=[1:T]
	I(t,:) = getI(Ipar,t*delta,nc);
end
noise = sigma*sqrt(delta)*randn(T,nc);

%-------------------------------------------------------------------------
% 			VOLTAGE TRACE
%-------------------------------------------------------------------------
fprintf('................ generating voltage trace \n')

% Use ode15s as solver. Only other one that might work well is ode23t, but twice as slow. 
sol =ode15s(@genkinetics,[delta Tmax],inivars',[],g,E,noise,delta,Tmax,C,R,W,nc,Ipar); 	

if usederiv==1;		[tg] = deval(sol,[delta:delta:Tmax]);
else			[tg,dv] = deval(sol,[delta:delta:Tmax]);
			dv=dv(1:nc,:)'*delta;
			dv=dv(:);
end
tg=tg';V=tg(:,1:nc);		% channel inference voltage

%................. derivative of voltage) .................................
if usederiv==1;
	for k=1:nc
		ind2 = (k-1)*T+1:k*T;
		dv(ind2,1) = gradient(V(:,k),delta)*delta;
	end
end

%................ macroscopic open fraction given voltage ................
for k=1:nc
	Jcgen(:,k     ) = (E(1)-V(:,k)) .* tg(:,k+nc   ).^3.*tg(:,k+2*nc );	
	Jcgen(:,k+  nc) = (E(2)-V(:,k)) .* tg(:,k+3*nc ).^4;			
	Jcgen(:,k+2*nc) = (E(3)-V(:,k));						
end

