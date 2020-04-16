%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of CONFIGR,
% as described in:
% Carpenter, G. A., Gaddam, C. S., & Mingolla, E. (2007). 
% CONFIGR: A vision-based system for long-range figure completion. 
% Neural Networks, xx(x) xxx-xxx. 
% Technical Report CAS/CNS TR-2007-016, Boston, MA: Boston University.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Chaitanya Sai (August 2007)
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%This function checks to see if an empty rectangle is bounded by tracks


function OutBool=CheckBound(LocImin,LocImax,LocJmin,LocJmax,SubPixRes,IcrawlT)

Horiz_=[LocJmin LocJmin+2:SubPixRes:LocJmax-2 LocJmax];
Vert_=[LocImin LocImin+2:SubPixRes:LocImax-2 LocImax];

Vert_tracks_A=sum(sign02(IcrawlT{2}(Vert_,LocJmin)+IcrawlT{4}(Vert_,LocJmin)));
Vert_tracks_B=sum(sign02(IcrawlT{2}(Vert_,LocJmax)+IcrawlT{4}(Vert_,LocJmax)));

Horiz_tracks_A=sum(sign02(IcrawlT{1}(LocImin,Horiz_)+IcrawlT{3}(LocImin,Horiz_)));
Horiz_tracks_B=sum(sign02(IcrawlT{1}(LocImax,Horiz_)+IcrawlT{3}(LocImax,Horiz_)));

if (Horiz_tracks_A+Horiz_tracks_B+Vert_tracks_A+Vert_tracks_B)==2*(size(Vert_,2)+size(Horiz_,2));
    OutBool=1;
else
    OutBool=0;
end

    