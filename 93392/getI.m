function I = getI(Ipar,t,nc);
%
% [I] = GETI(Ipar,t,nc);
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
% This function evaluates the squared sinusoidal current injected into 
% the cell at each point in time t. 
%
% Copyright Quentin Huys 2006


I = Ipar{1}+ Ipar{2}*sin(Ipar{3}*(repmat(t*2*pi,1,nc)+pi)).^2;
I(:,Ipar{4}) = 0;
