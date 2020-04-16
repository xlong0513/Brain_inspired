function A_Normalized = NormalizeA( A, normalizationMethod, aNorm)
% A_Normalized = NormalizeA( A,normalizationMethod,aNorm)
% This function normlizes the columns of A
% A: the matrix to be normalized
% normalizationMethod: 'L2 norm' or 'unit abs'; 
% aNorm: norm

% Set defaut value for 'normalizationMethod'
if ~exist('normalizationMethod', 'var')
    normalizationMethod = 'L2 norm';
end

% Set defaut value for 'aNorm'
if ~exist('aNorm', 'var')
    aNorm = 1;
end

% A small value to avoid zero division
epsilon = 1e-17;

% Normalize A
if isequal(normalizationMethod, 'L2 norm')
    % Normalize each column to l2-norm
    A_Normalized = A * diag( aNorm ./ ( sqrt( sum(A.*A,1) )+epsilon ) ); 
elseif isequal(normalizationMethod,'unit abs')   
    % Normalize each column such that the sum of absolute values of column
    % elements equals to 'aNorm'
    A_Normalized =  A * diag( aNorm ./ ( sum(abs(A),1) ) ); 
end

end