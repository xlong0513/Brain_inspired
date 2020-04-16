% MAIN.M
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
% This script is the main script. It calls PARAM.M to set the parameters of the
% cell to be simulated and used for inference. GETDATA.M simulates the cell and
% collects all the data needed to infer the parameters of interest. It calls
% GETI.M and also GENKINETCS.M as part of a call to a differential equation
% solver. SETUPQUAD.M sets up the quadratic programme to be solved. This is
% finally solved here.  PLOTS.M plots a small-scale version of figure 8 in the
% paper, and PLOTCELL.M finally plots the cells (you'll need some additional
% software for this, see in the file). 
%
% Copyright Quentin Huys 2006

clear all;
%........................ set up things ....................................
param;			% get parameters
getdata;		% get voltage trace, synaptic input, current etc. 
setupquad;		% set up the matrices M and the vector f

%........................ solve quadratic programme ..........................
fprintf('\n.............. solve quadratic programme\n')

switch minimizer
case 'qp'
	options = optimset('Display','none','MaxIter',5000);
	warning('off','optim:quadprog:SwitchToMedScale');
	warning('off','optim:quadprog:FullAndMedScale');

	a = quadprog(A,f,-eye(size(A)),zeros(size(f)));

case 'minq'
	a = minq(0,f,A,zeros(size(f)),repmat(Inf,size(f)),0);
end

plots;			% generate the plot in the paper

if doplotcell;
	plotcell;	
end
