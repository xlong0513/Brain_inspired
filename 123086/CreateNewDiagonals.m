
SubPixRes=5;
%load tempIpixelsMediumScaleNewRule
%load tempIpixelsCoarsestScaleNewRule
load MontereyScaledNew
%load MontereyScaledTwice

Iinterpol_Coarse=Iinterpol;
Ipixels_Coarse=Ipixels;
StoreDiagVals_Coarse=FillingStats.StoreDiagVals;
%filled_in_iter_Coars=filled_in_iter;

%load tempIpixelsFullScale

% load MontereyBig
% 
% IinterpolLines_Coarse=(zeros(size(Ipixels)));
% for i=1:size(StoreDiagVals_Coarse,1)
%     %CHANGE MULTIPLIER DEPENDING ON SCALE
%     tmpVect=[StoreDiagVals_Coarse(i,1:4)*2 StoreDiagVals_Coarse(i,5)];
%     IinterpolLines_Coarse(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))=InterpolateLines(IinterpolLines_Coarse(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4)),tmpVect(5),SubPixRes);
% end


IinterpolLines_Coarse=(zeros(size(Ipixels)*2));
for i=1:size(StoreDiagVals_Coarse,1)
    %CHANGE MULTIPLIER DEPENDING ON SCALE
    tmpVect=[StoreDiagVals_Coarse(i,1:4)*2 StoreDiagVals_Coarse(i,5)];
    IinterpolLines_Coarse(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4))=InterpolateLines(IinterpolLines_Coarse(tmpVect(1):tmpVect(2),tmpVect(3):tmpVect(4)),tmpVect(5),SubPixRes);
end
