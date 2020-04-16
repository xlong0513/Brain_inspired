function ImageColor=SimpleRoadImageBlack(Iin)

SizeIn=size(Iin);
Ibackground=(Iin==0);
Ifeature=(Iin==1);
IFfeature=(Iin==.5);%+(Iin==.75);
IFground=(Iin==.25)+(Iin==.75);
IEdge=(Iin==.6);
Igrids=[Iin==.375];
%Iconflict=(Iin==.75);
Idiagonal=(Iin==.85);
Icorner=(Iin==.125);
Ihighlight=(Iin==.542);
%This is the new version where the filled-in ground is darker
Grey2Val=[156+40 154+40 156+40]/255;
%GreyVal=[156+70 154+70 156+70]/255;
GreyVal=[1 1 1];
%Grey2Val=[1 1 1];
CornerVal=[10 10 10]/255;

% %This is the old version where the filled-in ground is Lighter
% GreyVal=[156+40 154+40 156+40]/255;
% Grey2Val=[156+70 154+70 156+70]/255;


%This is the Reddish Brown I have been using all along
RedVal=[174 69+.40 24+.50+.40]/255;


%This is the Black that Gail requested on July 21st
RedVal=[0.01 0.01 0.01];
HighlightVal=[255 140 0]/255;
%RedVal=[141+40 71+60 67+60]/255;
BrighterRedVal=[.75 0 0];
GreenVal=[222 186 18]/255;
DarkGreenVal=[74 146 8]/255;

% %This is actually Blue, used for coarser scale
GreenVal=[166 87 168]/255;
DarkGreenVal=[15 189 240]/255;

%DarkGreenVal=[15 189 240]/255;
%GreyVal=[134 134 134]/255;
%GreyVal=[1 1 1];

ConflictOrangeVal=[1 .5 0];
WhiteVal=[1 1 1];


Itemp=ones(SizeIn(1),SizeIn(2),3);
ImageColorT=GreyVal(1)*Itemp;

ImageColor(:,:,1)=Ibackground.*GreyVal(1)...
           +Ifeature.*RedVal(1)...
           +IFfeature.*GreenVal(1)...
           +IFground.*Grey2Val(1)...
           +Idiagonal.*DarkGreenVal(1)...
           +IEdge.*BrighterRedVal(1)...
           +Igrids.*WhiteVal(1)...
           +Icorner.*CornerVal(1)...
           +Ihighlight.*HighlightVal(1);     
           %+Iconflict.*ConflictOrangeVal(1);
       
ImageColor(:,:,2)=Ibackground.*GreyVal(2)...
           +Ifeature.*RedVal(2)...
           +IFfeature.*GreenVal(2)...
           +IFground.*Grey2Val(2)...
           +Idiagonal.*DarkGreenVal(2)...
           +IEdge.*BrighterRedVal(2)...
           +Igrids.*WhiteVal(2)...
           +Icorner.*CornerVal(2)...
           +Ihighlight.*HighlightVal(2);
           %+Iconflict.*ConflictOrangeVal(2);
           
ImageColor(:,:,3)=Ibackground.*GreyVal(3)...
           +Ifeature.*RedVal(3)...
           +IFfeature.*GreenVal(3)...
           +IFground.*Grey2Val(3)...
           +Idiagonal.*DarkGreenVal(3)...
           +IEdge.*BrighterRedVal(3)...
           +Igrids.*WhiteVal(3)...
           +Icorner.*CornerVal(3)...
           +Ihighlight.*HighlightVal(3);
           %+Iconflict.*ConflictOrangeVal(3);

 ImageColorT=repmat((ImageColor(:,:,1)==0).*(ImageColor(:,:,2)==0).*(ImageColor(:,:,3)==0),[1,1,3])*GreyVal(1);
 ImageColor=ImageColor+ImageColorT;
       
imagesc(ImageColor);
%SquareErase;