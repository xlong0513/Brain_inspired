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



tic




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
%Image files are stored in the following .mat files.
%%%%%%%%%

load CONFIGRtestimages;


%%%%%%%%%
%The Benchmark Monterey image is stored in the following .mat file
%%%%%%%%%
%load MontereyImage


%%%%%%%%%
%The default setting in the code is to run CONFIGR on the Monterey image
%%%%%%%%%

I=Ithick40;


%%%%%%%%%
%Set number of iterations and/or filling-limit here
%%%%%%%%%

ImSize=size(I)*5;
SubPixRes=5;
NumIter=round(SubPixRes*length(I)/2)+1; %% Change to 6 for size 201
FillNoFillVal=.25;
FillInLimit=round((250/max(ImSize))*max(ImSize));
FIGR_VAL=.5;
stp_diag_=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



I=imresize(I,ImSize);



%===These Lines are for size three subpixels========

I1=zeros(ImSize);
for i=1:ceil(ImSize(1)/SubPixRes)
    I1((i-1)*SubPixRes+1,:)=1;
end
for i=1:ceil(ImSize(2)/SubPixRes)
    I1(:,(i-1)*SubPixRes+1)=1;
end
% Itemp=zeros(201);
% Itemp(1:200,1:200)=I;
% I=Itemp; 

% ================================================== 



% Pixels are dilated to align them on the boundary grid. This ensures that
% boundaries can only occur at pre-ordained locations. The pixel dilation
% is necessary to ensure that figure pixels and boundary pixels are
% different. While in reality figure pixels are as thick as the resolution
% of the imaging device allows and boundary pixels infinitely thin, this
% cannot be achieved when manipulating images digitally

I=sign02((I+I1.*imdilate(I-I1.*I,[1 1 1;1 1 1;1 1 1])));
Ipixels=I;
IpixelsFillA=I;
IpixelsFillA_temp=sparse(I);
IpixelsFillA_real=sparse(I);

Iemergent_corners{1}=sparse(zeros(size(I)));
Iemergent_corners{2}=sparse(zeros(size(I)));
Iemergent_corners{3}=sparse(zeros(size(I)));
Iemergent_corners{4}=sparse(zeros(size(I)));

IpixelsTemp=Ipixels;
IpixelsFill=0*Ipixels;
IpixelsTempIter12=Ipixels;
IpixelsFillTempIter12=0*Ipixels;
IpixelsTempIter24=Ipixels;
IpixelsFillTempIter24=0*Ipixels;

IpixelsFilledInTemp=sparse(zeros(ImSize));
IpixelsFilledIn=sparse(zeros(ImSize));

Iinterpol=zeros(ImSize);
Icrawl=cell(5,1);
IgroundInters=ones(size(I))*1000;

IfillIter=zeros(size(I));

ItempCorners=zeros(ImSize);

IrectsShow=zeros(size(I));


rectsAll=1;
rectsAllNew=1;
resortArray=0