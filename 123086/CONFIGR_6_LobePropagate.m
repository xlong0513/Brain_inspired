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


%Initialize Temporary Matrices here
IpixelsFilledInTemp=sparse(zeros(ImSize));
Lob0=fliplr(LobeTemp1);
Lob180=(LobeTemp1);
Lob90=(LobeTemp1');
Lob270=flipud(LobeTemp1');




if kk~=2
    IcrawlT{1}=sign02((sign02(conv2(full(IcrawlT{1}),flipud(fliplr(Lob0)),'same'))... %Elongates any lobe activity at a corner in the image
        -sign02(conv2(full(IcrawlT{5})+full(IcrawlT{8}),flipud(fliplr(Lob0)),'same'))...% Newly formed corners stop propagation of contributing lobe activities
        ).*(1+Iedges-I) ...%Sides of pixels do-not contribute to lobe activity, Walls Block Propagation 
        +IcrawlT{1}).*(1+Iedges-I);
    
    IcrawlT{2}=sign02((sign02(conv2(full(IcrawlT{2}),flipud(fliplr(Lob90)),'same'))...
        -sign02(conv2(full(IcrawlT{5})+full(IcrawlT{6}),flipud(fliplr(Lob90)),'same'))...
        ).*(1+Iedges-I) ...
        +IcrawlT{2}).*(1+Iedges-I);
    
    IcrawlT{3}=sign02((sign02(conv2(full(IcrawlT{3}),flipud(fliplr(Lob180)),'same'))...
        -sign02(conv2(full(IcrawlT{6})+full(IcrawlT{7}),flipud(fliplr(Lob180)),'same'))...
        ).*(1+Iedges-I) ...
        +IcrawlT{3}).*(1+Iedges-I); 
    
    IcrawlT{4}=sign02((sign02(conv2(full(IcrawlT{4}),flipud(fliplr(Lob270)),'same'))...
        -sign02(conv2(full(IcrawlT{7})+full(IcrawlT{8}),flipud(fliplr(Lob270)),'same'))...
        ).*(1+Iedges-I) ...
        +IcrawlT{4}).*(1+Iedges-I);
else
    IcrawlT{1}=sign02((sign02(conv2(full(IcrawlT{1}),flipud(fliplr(Lob0)),'same'))... %Elongates any lobe activity at a corner in the image
        ).*(1+Iedges-I) ...%Sides of pixels do-not contribute to lobe activity, Walls Block Propagation 
        +IcrawlT{1}).*(1+Iedges-I);
    
    IcrawlT{2}=sign02((sign02(conv2(full(IcrawlT{2}),flipud(fliplr(Lob90)),'same'))...
        ).*(1+Iedges-I) ...
        +IcrawlT{2}).*(1+Iedges-I);
    
    IcrawlT{3}=sign02((sign02(conv2(full(IcrawlT{3}),flipud(fliplr(Lob180)),'same'))...
        ).*(1+Iedges-I) ...
        +IcrawlT{3}).*(1+Iedges-I); 
    
    IcrawlT{4}=sign02((sign02(conv2(full(IcrawlT{4}),flipud(fliplr(Lob270)),'same'))...
        ).*(1+Iedges-I)...
        +IcrawlT{4}).*(1+Iedges-I);
end


IcrawlTCheck=(IcrawlT{1}+IcrawlT{2}+IcrawlT{3}+IcrawlT{4});


IcrawlT{5}=sparse(sign02(IcrawlT{1}+IcrawlT{2}-1.1));
IcrawlT{6}=sign02(IcrawlT{2}+IcrawlT{3}-1.1);
IcrawlT{7}=sign02(IcrawlT{3}+IcrawlT{4}-1.1);
IcrawlT{8}=sign02(IcrawlT{4}+IcrawlT{1}-1.1);


if kk==1
    IshiftSW=conv2(I,[0 0 0;0 0 0;1 0 0],'same');
    IshiftSE=conv2(I,[0 0 0;0 0 0;0 0 1],'same');
    IshiftNE=conv2(I,[0 0 1;0 0 0;0 0 0],'same');
    IshiftNW=conv2(I,[1 0 0;0 0 0;0 0 0],'same');
    IcrawlT_img{1}=sparse(IcrawlT{1});
    IcrawlT_img{2}=sparse(IcrawlT{2});
    IcrawlT_img{3}=sparse(IcrawlT{3});
    IcrawlT_img{4}=sparse(IcrawlT{4});
    IcrawlT_img{5}=sparse(sign02(IcrawlT{1}+IcrawlT{2}-1.1).*IshiftSW);
    IcrawlT_img{6}=sparse(sign02(IcrawlT{2}+IcrawlT{3}-1.1).*IshiftSE);
    IcrawlT_img{7}=sparse(sign02(IcrawlT{3}+IcrawlT{4}-1.1).*IshiftNE);
    IcrawlT_img{8}=sparse(sign02(IcrawlT{4}+IcrawlT{1}-1.1).*IshiftNW);
end

IpixelsNon=sparse((Ipixels>0));
IshiftSW=conv2(double(full(IpixelsNon)),[0 0 0;0 0 0;1 0 0],'same');
IshiftSE=conv2(double(full(IpixelsNon)),[0 0 0;0 0 0;0 0 1],'same');
IshiftNE=conv2(double(full(IpixelsNon)),[0 0 1;0 0 0;0 0 0],'same');
IshiftNW=conv2(double(full(IpixelsNon)),[1 0 0;0 0 0;0 0 0],'same');


%EDIT: The addition of the shifted Ipixel helps incorporate the new
%empty and filled corner rules.

IcrawlTPura_temp{5}=sparse(sign02(IcrawlT{1}+IcrawlT{2}-1.1).*(1-IcrawlT_img{5})).*conv2(1-sign02(Ipixels),[0 0 0; 0 0 0; 1 0 0],'same');
IcrawlTPura_temp{6}=sparse(sign02(IcrawlT{2}+IcrawlT{3}-1.1).*(1-IcrawlT_img{6})).*conv2(1-sign02(Ipixels),[0 0 0; 0 0 0; 0 0 1],'same');
IcrawlTPura_temp{7}=sparse(sign02(IcrawlT{3}+IcrawlT{4}-1.1).*(1-IcrawlT_img{7})).*conv2(1-sign02(Ipixels),[0 0 1; 0 0 0; 0 0 0],'same');
IcrawlTPura_temp{8}=sparse(sign02(IcrawlT{4}+IcrawlT{1}-1.1).*(1-IcrawlT_img{8})).*conv2(1-sign02(Ipixels),[1 0 0; 0 0 0; 0 0 0],'same');

IcrawlTPura{5}=sparse(IcrawlTPura_temp{5}.*(1-sign02(IcrawlTPura_temp{6}+IcrawlTPura_temp{8})).*(1-IshiftSW));
IcrawlTPura{6}=sparse(IcrawlTPura_temp{6}.*(1-sign02(IcrawlTPura_temp{5}+IcrawlTPura_temp{7})).*(1-IshiftSE));
IcrawlTPura{7}=sparse(IcrawlTPura_temp{7}.*(1-sign02(IcrawlTPura_temp{6}+IcrawlTPura_temp{8})).*(1-IshiftNE));
IcrawlTPura{8}=sparse(IcrawlTPura_temp{8}.*(1-sign02(IcrawlTPura_temp{7}+IcrawlTPura_temp{5})).*(1-IshiftNW));


IcrawlTWall{5}=sparse(IcrawlTPura_temp{5}.*(sign02(IcrawlTPura_temp{6}+IcrawlTPura_temp{8})).*(1-IshiftSW));
IcrawlTWall{6}=sparse(IcrawlTPura_temp{6}.*(sign02(IcrawlTPura_temp{5}+IcrawlTPura_temp{7})).*(1-IshiftSE));
IcrawlTWall{7}=sparse(IcrawlTPura_temp{7}.*(sign02(IcrawlTPura_temp{6}+IcrawlTPura_temp{8})).*(1-IshiftNE));
IcrawlTWall{8}=sparse(IcrawlTPura_temp{8}.*(sign02(IcrawlTPura_temp{7}+IcrawlTPura_temp{5})).*(1-IshiftNW));
IcrawlAllWalls=sparse(sign02(IcrawlTWall{5}+IcrawlTWall{6}+IcrawlTWall{7}+IcrawlTWall{8}));


%This matrix stores 1's wherever emergent corners are created.
Iemergent_corners{1}=sparse(sign(Iemergent_corners{1}+IcrawlTPura{5}));
Iemergent_corners{2}=sparse(sign(Iemergent_corners{2}+IcrawlTPura{6}));
Iemergent_corners{3}=sparse(sign(Iemergent_corners{3}+IcrawlTPura{7}));
Iemergent_corners{4}=sparse(sign(Iemergent_corners{4}+IcrawlTPura{8}));




%Create Sparse Matrix: This was initially done to prevent the matrices from bloating unnecessarily as the matrix-size increased    
for i=1:4
    Isparse{i}=sparse(IcrawlT{4+i});
    [iSP{i},jSP{i},sSP{i}]=find(Isparse{i});
end

%sSPAll stores information of where the walls are located
for i=1:4
    IsparseAll{i}=sparse(IcrawlTPura_temp{4+i}+IcrawlTPura{4+i});
    [iSPAll{i},jSPAll{i},sSPAll{i}]=find(IsparseAll{i});
end


for i=1:4
    IsparsePura{i}=sparse(IcrawlTPura{4+i});
    [iSPPura{i},jSPPura{i},sSPPura{i}]=find(IsparsePura{i});
end

for i=1:4
    IsparseWall{i}=sparse(IcrawlTWall{4+i});
    [iSPWall{i},jSPWall{i},sSPWall{i}]=find(IsparseWall{i});
end


if kk==3*5+1
    
    'List of Empty FIGURE and GROUND Corners'
    [iSPPura{1},jSPPura{1},sSPPura{1}]
    
    [iSPPura{3},jSPPura{3},sSPPura{3}]
    
    
end
    
