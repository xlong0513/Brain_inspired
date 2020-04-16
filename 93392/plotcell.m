% PLOTCELL.M
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
% This script generates a plot of the cell if we are on a unix system
% this uses the GraphViz2Mat1.2 package, available from
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4518&objectType=file
%
% Copyright Quentin Huys 2006

if isunix 
	for k=1:nc; Lbl{k}=' ';end
	graph_to_dot(W,'directed',0,'node_label',Lbl);
	!neato -Tps -Gsplines=1 -Nheight=.01 -Nwidth=.01 -Nshape=box tmp.dot | cat > graph.ps ; gv graph.ps;
else error('This is not a UNIX system. Not sure this will work')
end

