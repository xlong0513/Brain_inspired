function y = makecolumn(x)

% makecolumn(x) returns a columnvector with the entries of x, if x is 
% vector. It returns an error message, if x i not a vector.
%
% created: Jakob Heinzle 01/07

 s = size(x);
 if min(s) > 1
     error('Input needs to be a vector');
 end
 if s(1)<s(2)
     y = x';
 else
    y = x;
 end
