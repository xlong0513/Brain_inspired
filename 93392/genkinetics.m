function [dV] = genkinetics(t,V,g,E,noise,delta,Tmax,C,R,W,nc,Ipar);
%
% [dV] = GENKINETICS.M
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
% This function simulates the cell under scrutiny. It is called within the
% script GETDATA.M. 
%
% Copyright Quentin Huys 2006


dV = zeros(nc*3,1);
I = getI(Ipar,t,nc);

%......................... overall equation..................................
for k=1:nc
dV(k ) = 1/C*( 	- g(k,1)*V(k+  nc)^3*V(k+2*nc)	*(V(k)-E(1)) ...
		- g(k,2)*V(k+3*nc)^4		*(V(k)-E(2)) ...
		- g(k,3)			*(V(k)-E(3)) ...
		+ R*I(k) + noise(ceil(t/delta),k) ...
		+ W(k,:)*(V(1:nc)-V(k)));
end

for k=1:nc

	%........................ gate constants ....................................
	%am1	=0.1*(V(k)+40)./(1-exp(-(V(k)+40)/10));		% Hodgkin Huxley gates
	%ah1	=0.07*exp(-(V(k)+65)/20);
	am1	=0.1*(V(k)+35)./(1-exp(-(V(k)+35)/10));		% Hodgkin Huxley gates -- modified
	ah1	=0.07*exp(-(V(k)+50)/20);
	an1	=0.01*(V(k)+55)./(1-exp(-(V(k)+55)/10));     
	bm1	=4*exp(-(V(k)+65)/18);
	bh1	=1./(exp(-(V(k)+35)/10)+1);
	bn1	=0.125*exp(-(V(k)+65)/80);

	%.......................... gate diffeq's.......................................

	dV(  nc+k)  = am1.*(1-V(  nc+k)) - bm1.*V(  nc+k);
	dV(2*nc+k)  = ah1.*(1-V(2*nc+k)) - bh1.*V(2*nc+k);
	dV(3*nc+k)  = an1.*(1-V(3*nc+k)) - bn1.*V(3*nc+k);
end

