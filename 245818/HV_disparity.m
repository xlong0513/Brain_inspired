function D = HV_disparity(II,n_scales,n_orient,energy_thr,th2,Dori)

%%% v0.1
%%% 23/10/2018
%%%
%%% M.Chessa and F. Solari
%%% University of Genoa, ITALY
%%%
%%% manuela.chessa@unige.it
%%% fabio.solari@unige.it
%%%
%%% REF PAPER:
%%% M. Chessa,  F. Solari. 
%%% A Computational Model for the Neural Representation and Estimation of the 
%%% Binocular Vector Disparity from Convergent Stereo Image Pairs. 
%%% International Journal of Neural Systems, 28, art. no. 1850029, 2018
%%% DOI: https://doi.org/10.1142/S0129065718500296
%%%
%%% INPUTS:
%%% II: image stereo pairs (:,:,1 is left ans  :,:,2 is right)
%%% n_scales: number of spatial scales
%%% n_orient: spatial orientations of V1 filters
%%% energy_thr: V1 global energy threshold 
%%% th2: V1  energy threshold 
%%% Dori: MT (disparity) directions
%%%
%%% OUTPUTS:
%%% D: computed horizontal and vertical disparity [m X n X 2]


suboct=1.414;%half-octave
II = image_pyramid(II,2,n_scales,suboct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimb=2;
ker=fspecial('average',dimb*2+1);
for sk=1:n_scales
    I=II{sk};
    [rr,cc,ff]=size(I);
    Im=zeros(rr,cc,ff);
    for kk=1:ff
        imb=putborde_img(I(:,:,kk),dimb,dimb);
        %mean
        seqm=conv2(imb,ker,'same');
        seqm(1:dimb,:)=[];
        seqm((end-dimb+1):end,:)=[];
        seqm(:,1:dimb)=[];
        seqm(:,(end-dimb+1):end)=[];
        %std
        seqs=stdfilt(imb,ones(5));
        seqs(1:dimb,:)=[];
        seqs((end-dimb+1):end,:)=[];
        seqs(:,1:dimb)=[];
        seqs(:,(end-dimb+1):end)=[];
        Im(:,:,kk)=(I(:,:,kk)-seqm)./(seqs+1e-10);
    end
    IIm{sk}=Im;
    %IIm{sk}=II{sk};%no (I-Im)/std
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpw=0.76;%disparity range
disp_values=-tmpw:(2*tmpw/(11-1)):tmpw;
fo=1/3.8;
for i=0:n_orient-1
    %theta=pi/2-(i*pi/8);
    ph_sh(i+1,:)=-2*pi*(fo)*disp_values;
end

%%%%%%%%%%%
% Level 1 %
%%%%%%%%%%%
F = filt_gabor_space(IIm{1},n_orient);

G{1,1}=F{1}(:,:,:,1);
G{1,2}=F{2}(:,:,:,1);
[G{2,1},G{2,2}]=shift_in_phase(F{1}(:,:,:,2),F{2}(:,:,:,2),ph_sh);
clear F
D = population(G,energy_thr,th2,disp_values,II{1},Dori);

invalid = isnan(D);
clear E
clear G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coarse-to-fine Estimation and Merging %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for scale = 2:n_scales
    
    [rr cc st] = size(II{scale});
    D = expand(D.*suboct,rr,cc);
    
    F = filt_gabor_space(IIm{scale},n_orient);
    G= warp_sequence(F,D);
    [G{2,1},G{2,2}]=shift_in_phase(G{2,1},G{2,2},ph_sh);
    
    clear F
    Ds = population(G,energy_thr,th2,disp_values,II{scale},Dori);
    
    clear G
    clear E
    D = merge_disparity(D,Ds);
    invalid = isnan(Ds);
    
end

D(invalid) = NaN;
%output
D=-D;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [II,x_pix,y_pix] = image_pyramid(II,n_frames,n_scales,suboct)


[sy, sx, nframe] = size(II);

tmp = II;
II = cell(1,n_scales);
II{n_scales} = tmp;

for scale = n_scales-1:-1:1
    
    for frame=1:nframe
        tmp=II{scale+1}(:,:,frame);
        tmp=imresize(tmp,1/suboct,'bilinear');
        II{scale}(:,:,frame) = tmp;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Oo = expand(O,rr,cc)

Oo(:,:,1) = imresize(O(:,:,1), [rr cc],'bilinear');
Oo(:,:,2) = imresize(O(:,:,2), [rr cc],'bilinear');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I2 = bilin_interp(I1,X,Y)
% function I2 = bilin_interp(I1,X,Y)
% Arbitrary image rewarping using bilinear interpolation
% X (column) and Y (row) (both floating point) are the source locations,
% used to fill the respective pixels

[nY1,nX1,rem] = size(I1); % source size
[nY2,nX2] = size(X); % target size

s = size(I1);
s(1:2) = [];

I2 = NaN.*zeros([ nY2 nX2 s ]);

for r = 1:rem
    
    for x = 1:nX2
        for y = 1:nY2
            
            % Pixel warping (2x2 group)
            
            x_w = floor(X(y,x));
            y_w = floor(Y(y,x));
            
            % Check validity
            
            if ( (x_w>0) && (x_w<nX1) && (y_w>0) && (y_w<nY1) )
                
                xs = X(y,x) - x_w;
                min_xs = 1-xs;
                ys = Y(y,x) - y_w;
                min_ys = 1-ys;
                
                w_00 = min_xs*min_ys;  % w_xy
                w_10 = xs*min_ys;
                w_01 = min_xs*ys;
                w_11 = xs*ys;
                
                I2(y,x,r) = w_00*I1(y_w,x_w,r) + w_10*I1(y_w,x_w+1,r) + ...
                    w_01*I1(y_w+1,x_w,r) + w_11*I1(y_w+1,x_w+1,r);
                
            end
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G_warped]= warp_sequence(G,D)
n_frame=2;
[sy sx n_orient phase_num re_im] = size(G{2});
G_warped = cell(2,n_frame);
[X Y] = meshgrid(1:sx,1:sy);

D(isnan(D)) = 0;


% keep first (left) frame

G_warped{1,1}=G{1,1}(:,:,:,1);
G_warped{1,2}=G{1,2}(:,:,:,1);

% warp second frame (right) to first (left)

Xn = X-D(:,:,1);
Yn = Y-D(:,:,2);
G_warped{2,1} = ...
    bilin_interp(double(G{1,1}(:,:,:,2)), Xn, Yn);
G_warped{2,2} = ...
    bilin_interp(double(G{1,2}(:,:,:,2)), Xn, Yn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function D = merge_disparity(D1,D2)


invalid1 = isnan(D1);
invalid2 = isnan(D2);

D1(invalid1) = 0;
D2(invalid2) = 0;

invalid = invalid1 & invalid2;

D = D1 + D2;
D(invalid) = NaN;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = population(G,th,th2,disp_values,II,Dori)


[sy sx n_orient n_phases re_im] = size(G{2});

CL = G{1,1};
CR = G{2,1};
SL = G{1,2};
SR = G{2,2};

E=zeros(sy,sx,n_orient,n_phases);

for o = 1:n_orient
    for ph_n=1:n_phases
        E(:,:,o,ph_n) = sqrt((CL(:,:,o) + squeeze(CR(:,:,o,ph_n))).^2 + (SL(:,:,o) + squeeze(SR(:,:,o,ph_n))).^2);
    end
end

E=E.^0.18;%V1 non-linearity
%
E=E/max(max(max(max(E))));
mask=E>th;
E=E.*mask;

D = NaN(sy, sx, 2);

%%%normalization: tmp matrix used later
tmp=zeros(sy,sx,n_phases);
for v=1:n_phases
    for o=1:n_orient
        tmp(:,:,v)=tmp(:,:,v)+E(:,:,o,v);
    end
end

tmpwg=11;%spatial support of pooling
filt = fspecial('gaussian', [tmpwg tmpwg], tmpwg/6);
filt=filt./sum(sum(fspecial('gaussian', [tmpwg tmpwg], tmpwg/6)));


for v=1:n_phases
    for o=1:n_orient
        E_v1(:,:,o,v)=E(:,:,o,v)./(tmp(:,:,v)+1e-9);%normalization
        E_v1(:,:,o,v) = conv2(E_v1(:,:,o,v), filt, 'same');%V1 spatial pooling
        %
        mask=ones(sy,sx);
        mask(tmp(:,:,v)<th2)=NaN;%threshold (unreliable pixels)
        E_v1(:,:,o,v)=E_v1(:,:,o,v).*mask;
    end
end

E_v1=0.81*E_v1/max(max(max(max(E_v1))));%exponential gain

%%%%%%orientation pooling
if Dori==2
    vdir=[0,pi/2];
else
    vdir=-5*pi/16:(10*pi/16)/(Dori-1):5*pi/16;
    vdir=vdir+pi/2;
end

for ivd=1:Dori
    for v=1:n_phases
        csum=zeros(sy,sx);
        for o=1:n_orient
            theta=(o-1)*(2*pi/n_orient);
            csum=csum + cos(vdir(ivd)-theta)*E_v1(:,:,o,v);
        end
        E_mt(:,:,ivd,v)=exp(csum);
    end
end

VD=zeros(sy,sx,Dori);
vx=zeros(sy,sx);
vy=zeros(sy,sx);


%%%%%%%%%interpolation: filling-in of borders and unreliable pixels
MMi=max(max(max(II(:,:,2))));
mmi=min(min(min(II(:,:,2))));

[xx1,xx2,xx3]=size(E_mt(:,:,1));

if xx1<20 && (~isequal(isnan(E_mt),zeros(size(E_mt))) )
    for ivd=1:Dori
        for v=1:n_phases
            Os=E_mt(:,:,ivd,v);
            Os=fillin_ppp(Os,II(:,:,2),(MMi-mmi)/3);
            E_mt(:,:,ivd,v)=Os;
        end
    end
end

if xx1>=20
    for ivd=1:Dori
        for v=1:n_phases
            Os=E_mt(:,:,ivd,v);
            Os(1:5,:,:)=NaN; Os((end-4):end,:,:)=NaN;  Os(:,1:5,:)=NaN; Os(:,(end-4):end,:)=NaN;
            Os=fillin_ppp(Os,II(:,:,2),(MMi-mmi)/3);
            E_mt(:,:,ivd,v)=Os;
        end
    end
end


%%%%%%MT decoding
%%%%component
for ivd=1:Dori
    for v=1:n_phases
        VD(:,:,ivd)=VD(:,:,ivd)+E_mt(:,:,ivd,v)*disp_values(v);
    end
end

%%%%argmin
for ivd=1:Dori
    vx=vx + VD(:,:,ivd)*cos(vdir(ivd));
    vy=vy + VD(:,:,ivd)*sin(vdir(ivd));
end
vx=(2/Dori)*vx;
vy=(2/Dori)*vy;

D(:,:,1)=vx; 
D(:,:,2)=vy; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IF = filt_gabor_space(I,n_filters)

[sy,sx,n_frames] = size(I);

IF{1} = zeros(sy,sx,n_filters,n_frames);
IF{2} = zeros(sy,sx,n_filters,n_frames);


w=5;%spatial support
f0=1/3.8;%peak frequency
sigma_s=1.2*sqrt(2*log(2)) ./ (2*pi*(f0/3));


[X,Y] = meshgrid(-w:w,-w:w);
theta=0:2*pi/n_filters:(2*pi-2*pi/n_filters);

G = exp(-(X.^2+Y.^2)/(2*sigma_s^2));

for ii=1:length(theta)
    XT=cos(theta(ii))*X+sin(theta(ii))*Y;
    GC=G.*cos(2*pi*f0*XT);
    GCB{ii}=GC-sum(sum(GC))/(2*w+1)^2;%DC
    GS=G.*sin(2*pi*f0*XT);
    GSB{ii}=GS-sum(sum(GS))/(2*w+1)^2;%DC
end


for frame = 1:n_frames
    
    for ii=1:n_filters/2
        even=conv2b(I(:,:,frame),GCB{ii});
        odd=conv2b(I(:,:,frame),GSB{ii});
        
        IF{1}(:,:,ii,frame) = even;
        IF{1}(:,:,ii+n_filters/2,frame) = even;
        
        IF{2}(:,:,ii,frame) = odd;
        IF{2}(:,:,ii+n_filters/2,frame) = -odd;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function imf=conv2b(im, ker)
[nky,nkx]=size(ker);
sh='valid';
Bx=(nkx-1)/2;
By=(nky-1)/2;
im=putborde(im,Bx,By);
imf=conv2(im,ker,sh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function imb=putborde(im,Nx,Ny)

[sy,sx]=size(im);
imb=zeros(sy+2*Ny,sx+2*Nx);
imb(1+Ny:sy+Ny,1+Nx:sx+Nx)=im;

for k=1:Nx
    imb(Ny+1:sy+Ny,k)=im(:,1);
    imb(Ny+1:sy+Ny,k+sx+Nx)=im(:,sx);
end
for k=1:Ny
    imb(k,Nx+1:sx+Nx)=im(1,:);
    imb(k+sy+Ny,Nx+1:sx+Nx)=im(sy,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OO=fillin_ppp(Oin,II,sigma_r)
%Boundary conditions&Unreliable regions

w=7;
sigma_d=(2*w+1)/6;
% Pre-compute Gaussian distance weights.
[X,Y] = meshgrid(-w:w,-w:w);
G = exp(-(X.^2+Y.^2)/(2*sigma_d^2));

tmp=Oin; tmp2=Oin;
tmp2(isnan(tmp2))=0;
masknonan=~isnan(tmp);
trueborder = bwmorph(masknonan,'remove');
tmp((trueborder)==0)=0;

% Apply bilateral filter.
dim = size(II);
tmpi=II;
%tmpi((trueborder)==0)=0;
B = zeros(dim);
for i = 1:dim(1)
    for j = 1:dim(2)
        
        % Extract local region.
        iMin = max(i-w,1);
        iMax = min(i+w,dim(1));
        jMin = max(j-w,1);
        jMax = min(j+w,dim(2));
        I = tmpi(iMin:iMax,jMin:jMax);
        O = tmp(iMin:iMax,jMin:jMax);
        
        % Compute Gaussian intensity weights.
        H = exp(-((I-tmpi(i,j)).^2)/(2*sigma_r^2));
        
        % Calculate bilateral filter response.
        GG=G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1);
        F = H.*GG;
        
        N=F.*trueborder(iMin:iMax,jMin:jMax);
        if sum(N(:))==0
            B(i,j) = sum(F(:).*O(:));
        else
            B(i,j) = sum(F(:).*O(:))/sum(N(:));
        end
        
    end
end

O1filled=(B).*(~masknonan) + tmp2.*masknonan;

OO=O1filled;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imb=putborde_img(im,Nx,Ny)

[sy,sx]=size(im);
imb=zeros(sy+2*Ny,sx+2*Nx);
imb(1+Ny:sy+Ny,1+Nx:sx+Nx)=im;

for k=1:Nx
    imb(Ny+1:sy+Ny,k)=im(:,1);
    imb(Ny+1:sy+Ny,k+sx+Nx)=im(:,sx);
end

for k=1:Ny
    imb(k,Nx+1:sx+Nx)=im(1,:);
    imb(k+sy+Ny,Nx+1:sx+Nx)=im(sy,:);
end

imb(1:Ny,1:Nx)=im(1,1);
imb(sy+Ny+1:end,1:Nx)=im(sy,1);
imb(1:Ny,sx+Nx+1:end)=im(1,sx);
imb(sy+Ny+1:end,sx+Nx+1:end)=im(sy,sx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, B]=shift_in_phase(F_real,F_imag,p)

n_orient=size(F_real,3);
[n_phase]=size(p);

for orient = 1:n_orient
    G_even_tmp=F_real(:,:,orient);
    G_odd_tmp=F_imag(:,:,orient);
    phase=1;
    
    for(ph=p(orient,:))
        G_even(:,:,orient,phase)=G_even_tmp*cos(ph)-G_odd_tmp*sin(ph);
        G_odd(:,:,orient,phase)=G_odd_tmp*cos(ph)+G_even_tmp*sin(ph);
        phase=phase+1;
    end
end
A = G_even;
B = G_odd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
