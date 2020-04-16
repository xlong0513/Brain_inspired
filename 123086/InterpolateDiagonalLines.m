%Usage:Iinterp = InterpolateDiagonalLines(I,Ipixels)
%This function interpolates the diagonals in CONFIGR filling-in areas
%I is the input to CONFIGR and Ipixels is the output obtained from
%[I_out,Ipixels]=runCONFIGR(I)


function Iinterp = InterpolateDiagonalLines(I,Ipixels)

CheckCornKernel=[-1 -1 -1;-1 1 1;-1 1 1];
CheckThinCornKernel=[-1 -1 -1;-1 1 -1;-1 1 -1];
IpixelsFigure=double((Ipixels==.5)+-1*(Ipixels~=.5));
IrectCorners=(conv2(IpixelsFigure,flipud(CheckCornKernel),'same')==9)+...
    (conv2(IpixelsFigure,fliplr(CheckCornKernel),'same')==9)+...
    (conv2(IpixelsFigure,flipud(fliplr(CheckCornKernel)),'same')==9)+...
    (conv2(IpixelsFigure,CheckCornKernel,'same')==9);
IrectCorners=IrectCorners+(conv2(IpixelsFigure,flipud(CheckThinCornKernel),'same')==9)+...
    (conv2(IpixelsFigure,rot90(CheckThinCornKernel),'same')==9)+...
    (conv2(IpixelsFigure,fliplr(rot90(CheckThinCornKernel)),'same')==9)+...
    (conv2(IpixelsFigure,CheckThinCornKernel,'same')==9);
IrectCornersT=double((conv2(I,[1 1 1;1 1 1; 1 1 1;],'same')==1).*IrectCorners);


[corn_i,corn_j]=find(IrectCornersT);

Iinterp=zeros(size(Ipixels));
for i=1:length(corn_i)
    for j=i:length(corn_i)
        if sum(sum(IpixelsFigure(min(corn_i(i),corn_i(j)):max(corn_i(i),corn_i(j)),min(corn_j(i),corn_j(j)):max(corn_j(i),corn_j(j)))~=1))==0
            try
            max_diff=max(abs(corn_i(i)-corn_i(j)),abs(corn_j(i)-corn_j(j)));
            if (corn_i(i)~=corn_i(j)) && (corn_j(i)~=corn_j(j))
            Iinterp(sub2ind(size(Iinterp),...
                   round(corn_i(i):sign(corn_i(j)-corn_i(i))*max((abs(corn_i(j)-corn_i(i)))/max_diff,.01):corn_i(j)),...
                   round(corn_j(i):sign(corn_j(j)-corn_j(i))*max((abs(corn_j(j)-corn_j(i)))/max_diff,.01):corn_j(j))))=1;
            else
                if abs(corn_i(i)-corn_i(j))>0
                    Iinterp(corn_i(i):sign(corn_i(j)-corn_i(i)):corn_i(j),corn_j(i))=1;
                else
                    Iinterp(corn_i(i),corn_j(i):sign(corn_j(j)-corn_j(i)):corn_j(j))=1;
                end
            end
            catch
                disp('some error')
            end
        end
    end
end
