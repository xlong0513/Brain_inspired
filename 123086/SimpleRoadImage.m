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


function ImageColor=SimpleRoadImage(Iin)

SizeIn=size(Iin);
Ibackground=(Iin==0);
Ifeature=(Iin==1);
IFfeature=(Iin==.5);
IFground=(Iin==.25)+(Iin==.75);
IEdge=(Iin==.6);
%Iconflict=(Iin==.75);
Idiagonal=(Iin==.85);
GreyVal=[.45 .45 .45];
RedVal=[.5 0 0];
BrighterRedVal=[.75 0 0];
GreenVal=[0.5 1 .5];
DarkGreenVal=[0 .75 0];
Grey2Val=[.6 .6 .6];
ConflictOrangeVal=[1 .5 0];


Itemp=ones(SizeIn(1),SizeIn(2),3);
ImageColorT=GreyVal(1)*Itemp;

ImageColor(:,:,1)=Ibackground.*GreyVal(1)...
           +Ifeature.*RedVal(1)...
           +IFfeature.*GreenVal(1)...
           +IFground.*Grey2Val(1)...
           +Idiagonal.*DarkGreenVal(1)...
           +IEdge.*BrighterRedVal(1);
           %+Iconflict.*ConflictOrangeVal(1);
       
ImageColor(:,:,2)=Ibackground.*GreyVal(2)...
           +Ifeature.*RedVal(2)...
           +IFfeature.*GreenVal(2)...
           +IFground.*Grey2Val(2)...
           +Idiagonal.*DarkGreenVal(2)...
           +IEdge.*BrighterRedVal(2);
           %+Iconflict.*ConflictOrangeVal(2);
           
ImageColor(:,:,3)=Ibackground.*GreyVal(3)...
           +Ifeature.*RedVal(3)...
           +IFfeature.*GreenVal(3)...
           +IFground.*Grey2Val(3)...
           +Idiagonal.*DarkGreenVal(3)...
           +IEdge.*BrighterRedVal(3);
           %+Iconflict.*ConflictOrangeVal(3);

 ImageColorT=repmat((ImageColor(:,:,1)==0).*(ImageColor(:,:,2)==0).*(ImageColor(:,:,3)==0),[1,1,3])*GreyVal(1);
 ImageColor=ImageColor+ImageColorT;
       
imagesc(ImageColor);
%SquareErase;