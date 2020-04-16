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



%This file contains the initialization commands; Set input location here.
CONFIGR_6_Init_Function;

%Initialize matrices
grnd_step=1;
figr_step=1;
rule3_step=1;
filled_in_grnd_rects_store=[0 0 0 0 0];
filled_in_rule3_rects_store=[0 0 0 0 0];
filled_in_figr_rects_store=[0 0 0 0 0];
fig_num_in=1;


%This file creates boundaries the first time; simple and complex cells
CONFIGR_6_FindBound;

%Lobe Propagation
for kk=1:NumIter
    %This file computes lobe propagation
    CONFIGR_6_LobePropagate;
    
    %This loop finds corner pairs for upper-right |_ and lower-left corners
    
    if mod(kk-1,5)==0  
        
        %This Ccunts the number of ADDITIONAL fillings-in due to Rule 3;
        NumFilled_Rule3(round((kk-1)/5)+1)=0;
        NumFilled_GROUND(round((kk-1)/5)+1)=0;
        NumFilled_FIGURE(round((kk-1)/5)+1)=0;
        
        
        
        EmptyRectangleTypeOne_Ground_6
        
        disp('Iteration Number')
        disp(round((kk-1)/5)+1)
        
        
        EmptyRectangleTypeTwo_Ground_6
        
        
        
        %Invoke code to star filling-in as GROUND
        FillingGROUND
        
        %Invoke code to star filling-in as GROUND
        FillingFIGURE
        
        %Can ignore: This is a porition of code involved in memory clean-up
        Inonsides=sign02(Inonsides-sign02(double((IcrawlT{1}+IcrawlT{3})>1)+double((IcrawlT{2}+IcrawlT{4})>1)));
        
        %Plot resultant CONFIGR output
        %plotIpixels(Ipixels,Iinterpol,filled_in_figr_rects_store,fig_num_in)
        fig_num_in=fig_num_in+1;
    end
    
end


%Store resultant fillings-in and corner locations for debugging
% 
% FillingStats.Rule=NumFilled_Rule3;
% FillingStats.GROUND=NumFilled_GROUND;
% FillingStats.FIGURE=NumFilled_FIGURE;
% FillingStats.rule3_rects=filled_in_rule3_rects_store;
% FillingStats.figr_rects=filled_in_figr_rects_store;
% FillingStats.grnd_rects=filled_in_grnd_rects_store;
% FillingStats.StoreDiagVals=StoreDiagVals;

%Use this to plot a snippet of the output image
%CONFIGR_6_CreatePlots;