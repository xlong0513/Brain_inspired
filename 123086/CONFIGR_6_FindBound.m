
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of CONFIGR,
% as described in CAS/CNS Technical Report TR-2007-016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Chaitanya Sai (August 2007)
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Finds boundaries in input image
for i=1:4
    Icrawl{i}=zeros(ImSize);
end


KernVert=[-1 .2 1];
KernHoriz=[-1 .2 1]';
Icrawl{1}=sparse(double(abs(conv2(I,KernHoriz,'same')-.2)==1));
Icrawl{2}=sparse(double(abs(conv2(I,KernVert,'same')-.2)==1));
Icrawl{3}=sparse(double(abs(conv2(I,KernHoriz,'same')-.2)==1));
Icrawl{4}=sparse(double(abs(conv2(I,KernVert,'same')-.2)==1));
Icrawl{9}=sparse(double(sign02(Icrawl{1}+Icrawl{2}+Icrawl{3}+Icrawl{4}-2)));
Inonsides=1+Icrawl{9}-I;

% Any non zero fractional value can be used instead of the .2 
KernVertA=[-1 .2 1;-1 .2 1;-1 .2 1];
KernVertB=[1 .2 -1;1 .2 -1;1 .2 -1];

KernHorizA=KernVertA';
KernHorizB=KernVertB';

KernCorn{1}=[-1 1 0];
KernCorn{2}=[0 1 -1]';
KernCorn{3}=[0 1 -1];
KernCorn{4}=[-1 1 0]';


Icrawl{1}=sparse(double(conv2(I,KernHorizA,'same')>2.5).*I+double(conv2(I,KernHorizB,'same')>2.5).*I);
Icrawl{2}=sparse(double(conv2(I,KernVertA,'same')>2.5).*I+double(conv2(I,KernVertB,'same')>2.5).*I);
Icrawl{3}=Icrawl{1};
Icrawl{4}=Icrawl{2};
% Icrawl{1}--0 degrees,Icrawl{2}--90,Icrawl{1}--180 degrees,Icrawl{1}--270 degrees


Icrawl{9}=sparse(double(sign02(Icrawl{1}+Icrawl{2}+Icrawl{3}+Icrawl{4}-2)));
Iedges=sign02(Icrawl{1}+Icrawl{2}+Icrawl{3}+Icrawl{4});
Iedges=I-imerode(I,[1 1 1;1 1 1;1 1 1]);

% These matrices store ones in the locations of the corners in the input
% image

Icrawl{5}=sparse(sign02(Icrawl{1}+Icrawl{2}-1.1)*(complex(1,1)));
Icrawl{6}=sparse(sign02(Icrawl{2}+Icrawl{3}-1.1)*(complex(-1,1)));
Icrawl{7}=sparse(sign02(Icrawl{3}+Icrawl{4}-1.1)*(complex(-1,-1)));
Icrawl{8}=sparse(sign02(Icrawl{4}+Icrawl{1}-1.1)*(complex(1,-1)));

% Inonsides stores one in the locations where there are no non-convex-corner sub-pixels.
%The non-convex specificity is a fortuitous outcome of the method used in computing boundary pixels    



for i=1:4;
    Icrawl{i}=Icrawl{i}-Icrawl{9}.*sign02(conv2(full(Icrawl{i}),KernCorn{i},'same'));
end

%Bipole Lobe Size
LobSiz=SubPixRes*2+1; 



LobeTemp1=zeros(LobSiz);
LobeTemp1(ceil(LobSiz/2),ceil(LobSiz/2)-1:ceil(LobSiz/2))=1;



IcrawlT{1}=sparse(Icrawl{1}.*(1+Iedges-I));
IcrawlT{2}=sparse(Icrawl{2}.*(1+Iedges-I));
IcrawlT{3}=sparse(Icrawl{3}.*(1+Iedges-I));
IcrawlT{4}=sparse(Icrawl{4}.*(1+Iedges-I));
IcrawlT{5}=sparse(zeros(size(IcrawlT{1})));
IcrawlT{6}=sparse(IcrawlT{5});
IcrawlT{7}=sparse(IcrawlT{5});
IcrawlT{8}=sparse(IcrawlT{5});

IcrawlT_old=IcrawlT;

IcrawlTPura{5}=sparse(zeros(size(IcrawlT{1})));
IcrawlTPura{6}=sparse(IcrawlT{5});
IcrawlTPura{7}=sparse(IcrawlT{5});
IcrawlTPura{8}=sparse(IcrawlT{5});

IcrawlTPura_temp{5}=sparse(zeros(size(IcrawlT{1})));
IcrawlTPura_temp{6}=sparse(IcrawlT{5});
IcrawlTPura_temp{7}=sparse(IcrawlT{5});
IcrawlTPura_temp{8}=sparse(IcrawlT{5});
