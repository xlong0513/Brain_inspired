function X = logistic(X)
% Logistic function.
%   Florian Raudies, 01/30/2014, Boston University.
X = 1./(1 + exp(-X));