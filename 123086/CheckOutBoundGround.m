%This function checks to see if a rectangle is bounded by tracks
function [OutBool , NumFilledOut]=CheckOutBoundGround(LocIminT,LocImaxT,LocJminT,LocJmaxT,SubPixRes,IpixelsFillA,FillNoFillVal,ImSize,NumFilled)


LocImin=max(LocIminT-1,1);
LocImax=min(LocImaxT,ImSize(1));

LocJmin=max(LocJminT-1,1);
LocJmax=min(LocJmaxT,ImSize(2));

Horiz_=[LocJminT+2:SubPixRes:LocJmaxT-2];
Vert_=[LocIminT+2:SubPixRes:LocImaxT-2];

%To remove Rule 3 (More than Half pixels are filled in or image FiGURE),
%simply make Track_length_half real big

Track_length_half=(length(Horiz_)+length(Vert_))*max(ImSize)*2;

%CHOOSE ONE OF THESE--|^

%Track_length_half=(length(Horiz_)+length(Vert_));

Vert_tracks_A=sum(IpixelsFillA(Vert_,LocJmin)==FillNoFillVal)+sum(IpixelsFillA(Vert_,LocJmax)==FillNoFillVal);
Horiz_tracks_A=sum(IpixelsFillA(LocImin,Horiz_)==FillNoFillVal)+sum(IpixelsFillA(LocImax,Horiz_)==FillNoFillVal);

Vert_tracks_Fig=sum(IpixelsFillA(Vert_,LocJmin)==.5)+sum(IpixelsFillA(Vert_,LocJmax)==.5)+sum(IpixelsFillA(Vert_,LocJmin)==1)+sum(IpixelsFillA(Vert_,LocJmax)==1);
Horiz_tracks_Fig=sum(IpixelsFillA(LocImin,Horiz_)==.5)+sum(IpixelsFillA(LocImax,Horiz_)==.5)+sum(IpixelsFillA(LocImin,Horiz_)==1)+sum(IpixelsFillA(LocImax,Horiz_)==1);

if (Horiz_tracks_A+Vert_tracks_A)>0 | ((Horiz_tracks_Fig+Vert_tracks_Fig)>Track_length_half);
    OutBool=1;
else
    OutBool=0;
end

NumFilledOut=NumFilled;
if (Horiz_tracks_A+Vert_tracks_A)==0 && ((Horiz_tracks_Fig+Vert_tracks_Fig)>Track_length_half);
NumFilledOut=NumFilled+1;

end
