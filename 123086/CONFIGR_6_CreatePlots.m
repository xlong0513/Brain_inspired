%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of CONFIGR,
% as described in:
% Carpenter, G. A., Gaddam, C. S., & Mingolla, E. (2007). 
% CONFIGR: A vision-based system for long-range figure completion. Neural Networks, xx(x) xxx-xxx. 
%  Technical Report CAS/CNS TR-2007-016, Boston, MA: Boston University.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Chaitanya Sai (August 2007)
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FillingStats.Rule=NumFilled_Rule3;
FillingStats.GROUND=NumFilled_GROUND;
FillingStats.FIGURE=NumFilled_FIGURE;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting CONFIGR outputs along with all the corners and lobe activities can be memory intensive. 
%This code extracts a part of the image centered on the location input
%below.
                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CenterPX=125;
CenterPY=125;
% 
% 
% 
Igrids=zeros(ImSize);
TempG_a=[1:SubPixRes:ImSize(2)];
TempG_b=[1:SubPixRes:ImSize(1)];
Igrids(:,TempG_a)=1;
Igrids(TempG_b,:)=1;



InpBipT{1}(:,:,1)=IcrawlT{1}(max(CenterPX-99,1):min(CenterPX+100,ImSize(2)),max(CenterPY-99,1):min(CenterPY+100,ImSize(2)));
InpBipT{1}(:,:,2)=IcrawlT{2}(max(CenterPX-99,1):min(CenterPX+100,ImSize(2)),max(CenterPY-99,1):min(CenterPY+100,ImSize(2)));
InpBipT{1}(:,:,3)=IcrawlT{3}(max(CenterPX-99,1):min(CenterPX+100,ImSize(2)),max(CenterPY-99,1):min(CenterPY+100,ImSize(2)));
InpBipT{1}(:,:,4)=IcrawlT{4}(max(CenterPX-99,1):min(CenterPX+100,ImSize(2)),max(CenterPY-99,1):min(CenterPY+100,ImSize(2)));
Iview=Ipixels(max(CenterPX-99,1):min(CenterPX+100,ImSize(2)),max(CenterPY-99,1):min(CenterPY+100,ImSize(2)));
IgridView=Igrids(max(CenterPX-99,1):min(CenterPX+100,ImSize(2)),max(CenterPY-99,1):min(CenterPY+100,ImSize(2)));



IinterpolV=Iinterpol(max(CenterPX-99,1):min(CenterPX+100,ImSize(2)),max(CenterPY-99,1):min(CenterPY+100,ImSize(2)));
IrectsShowV=IrectsShow(max(CenterPX-99,1):min(CenterPX+100,ImSize(2)),max(CenterPY-99,1):min(CenterPY+100,ImSize(2)));

toc
