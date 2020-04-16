function S = num2cellstr(N)
% num2cellstr
%   N       - Integer numbers.
%
% RETURN
%   S       - Cell array with numbers as strings.
%
%
%   Florian Raudies, 01/30/2014, Boston University.

nNum = numel(N);
S    = cell(nNum,1);
for iNum = 1:nNum, S{iNum} = num2str(N(iNum)); end