function A=sign01(B);
[sizBR sizBC]=size(B);
for i=1:sizBR
    for j=1:sizBC
A(i,j)=sign(B(i,j));
if A(i,j)==0;A(i,j)=1;end;
end
end


    