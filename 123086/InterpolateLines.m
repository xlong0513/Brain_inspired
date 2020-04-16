%Created on 11/14/2005
%This Function creates the interpolating diagonals of the empty rectangles
%that span figure in SimpleCrawlMayConflict. The algorithm is general
%enough that it creates a diagonal in any rectangle. type decides which
%main diagonal is to be spanned. type='normal' chooses the identity
%diagonal. 
%   Algorithm: 
%   1 Take Rectangle as Input; Compute Smallest size (row or column)
%   2 Compute an eye Matrix of that size and resize it size of Input.
%   3 Xor this with input and present as output.
%   Note: Input is assumed to be binary. Output is constrained to be binary
function Iout = InterpolateLines(Iin,type,SubPixRes)
IinSize=size(Iin);
minIinSize=min(IinSize);
Ieye=eye(round((minIinSize-1)/SubPixRes));
if type==1
Iout=sign02(imresize_old(Ieye,IinSize)+Iin);
else
    Iout=sign02(imresize_old(fliplr(Ieye),IinSize)+Iin);
end