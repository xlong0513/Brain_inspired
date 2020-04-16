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
%
%
% Usage:
% I_output = runCONFIGR(I,PixRes)
% I ---------- input image
% PixRes ----- pixel resolution 1:fine 2:medium 3:coarse
%              N implies that each coarse-scale pixel is created from an
%              2^Nx2^N square of fine-scale (original) pixels
% This runs CONFIGR with the following defaults:
% PixRes: pixel resolution default=1: CONFIGR pixel resolution is the same as that of the input image
%
%%%%%%%Options:%%%%%%%%%%%%%%%%
% I_output = runCONFIGR(I,PixRes,NumIter,ShrinkFact)
% I ----------- input image
% PixRes ------ pixel resolution 1:fine 2:medium 3:coarse
% NumIter ----- number of iterations CONFIGR simulation can be forced to stop early by setting a low number of iterations.
% ShrinkFact -- ratio of desired image size to actual size Sparse images can be resized for faster runtimes.
% Raw CONFIGR output (with ground and figure-filled rectangles)
% can be obtained using [I_output, I_output_raw, Idiagonals]=runCONFIGR(I)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_output,varargout] = runCONFIGR(I,varargin)


%Use double precision for I

I=double(I);

if ndims(I)~=2
    eval('help runCONFIGR')
    error('Error using CONFIGR code: Input must be a two-dimensional binary image')
elseif sum((I(:)~=0).*(I(:)~=1))~=0
    eval('help runCONFIGR')
    error('Error using CONFIGR code: Input must be a two-dimensional binary image')
end



if nargout>3
    eval('help runCONFIGR')
    error('Error using CONFIGR code: Only one,two, or three outputs allowed')
end

if (nargin<1) || (nargin>4)
    eval('help runCONFIGR')
    error('Error using CONFIGR code:  incorrect number of inputs. See directions above')
end


if nargin==1
    PixRes=1;
    NumIter=-1;
    ShrinkFact=1;
elseif nargin==2
    PixRes=round(varargin{1});

    if PixRes==0
        eval('help runCONFIGR')
        error('Error using CONFIGR: pixel resolution must be an integer greater than 0')
    end
    NumIter=-1;
    ShrinkFact=1;
elseif nargin==3
    PixRes=round(varargin{1});
    if PixRes==0
        eval('help runCONFIGR')
        error('Error using CONFIGR: pixel resolution must be an integer greater than 0')
    end
    NumIter=round(varargin{2})*5;
    if NumIter<0
        eval('help runCONFIGR')
        error('Error using CONFIGR: number of iterations must be an integer greater than 0')
    end
    ShrinkFact=1;
elseif nargin==4
    PixRes=round(varargin{1});
    if PixRes==0
        eval('help runCONFIGR')
        error('Error using CONFIGR: pixel resolution must be an integer greater than 0')
    end
    NumIter=round(varargin{2})*5;
    if NumIter<0
        eval('help runCONFIGR')
        error('Error using CONFIGR: number of iterations must be an integer greater than 0')
    end
    ShrinkFact=round(varargin{3});
end


if PixRes~=1
    for i=1:PixRes-1
        Itemp=double(conv2(I,[1/2 1/2 1/2 1/2;1/2 1 1 1/2; 1/2 1 1 1/2;1/2 1/2 1/2 1/2],'same')>=3);
        I=Itemp(1:2:size(Itemp,1),1:2:size(Itemp,2));
    end
end

if ShrinkFact~=1
    Itemp=double(conv2(I,ones(ShrinkFact),'same')>=1);
    I=Itemp(1:ShrinkFact:size(Itemp,1),1:ShrinkFact:size(Itemp,2));
end

if NumIter==-1
    NumIter=length(I)*5;
end


SubPixRes=5;
FillNoFillVal=.25;
FIGR_VAL=.5;
DIAG_VAL=.85;
ImSizeOrig=size(I);
ImSize=size(I)*5;
FillInLimit=round((250/max(ImSize))*max(ImSize));
stp_diag_=1;
ver_mat=ver;
ver_mat=str2num(ver_mat(strcmp('MATLAB',{ver_mat.Name})).Version(1:3))<7;




try

    try
        if ~ver_mat
            I=imresize_old_CONFIGR(I,ImSize);
        else
            I=imresize(I,ImSize);
        end
    catch
        error('The function imresize (Version 6.5) or imresize_old (Version 7 or greater) from the image toolbox is required')
    end



    CONFIGR_6

    if ~ver_mat
        I_output=imresize_old_CONFIGR(double(I+Iinterpol.*(1-I)),ImSizeOrig);

        if nargout==2
            varargout(1)={double(imresize_old_CONFIGR(double(Ipixels==1),ImSizeOrig)+...
                FIGR_VAL*imresize_old_CONFIGR(double(Ipixels==FIGR_VAL),ImSizeOrig)+ ...
                FillNoFillVal*imresize_old_CONFIGR(double(Ipixels==FillNoFillVal),ImSizeOrig)+ ...
                (DIAG_VAL-FIGR_VAL)*imresize_old_CONFIGR(Iinterpol.*(1-I),ImSizeOrig))};
        end
        if nargout==3
            varargout(1)={double(imresize_old_CONFIGR(double(Ipixels==1),ImSizeOrig)+...
                FIGR_VAL*imresize_old_CONFIGR(double(Ipixels==FIGR_VAL),ImSizeOrig)+ ...
                FillNoFillVal*imresize_old_CONFIGR(double(Ipixels==FillNoFillVal),ImSizeOrig)+ ...
                (DIAG_VAL-FIGR_VAL)*imresize_old_CONFIGR(Iinterpol.*(1-I),ImSizeOrig))};

            varargout(2)={double(imresize_old_CONFIGR(Iinterpol.*(1-I),ImSizeOrig))};
        end

    else
        I_output=imresize(double(I+Iinterpol.*(1-I)),ImSizeOrig);

        if nargout==2
            varargout(1)={double(imresize(double(Ipixels==1),ImSizeOrig)+...
                FIGR_VAL*imresize(double(Ipixels==FIGR_VAL),ImSizeOrig)+ ...
                FillNoFillVal*imresize(double(Ipixels==FillNoFillVal),ImSizeOrig)+ ...
                (DIAG_VAL-FIGR_VAL)*imresize(Iinterpol.*(1-I),ImSizeOrig))};
        end

        if nargout==3
            varargout(1)={double(imresize(double(Ipixels==1),ImSizeOrig)+...
                FIGR_VAL*imresize(double(Ipixels==FIGR_VAL),ImSizeOrig)+ ...
                FillNoFillVal*imresize(double(Ipixels==FillNoFillVal),ImSizeOrig)+ ...
                (DIAG_VAL-FIGR_VAL)*imresize(Iinterpol.*(1-I),ImSizeOrig))};

            varargout(2)={imresize(double(Iinterpol.*(1-I),ImSizeOrig))};
        end

    end

    disp('To plot output: SimpleRoadImage(I_output); To plot output with black input pixels: SimpleRoadImageBlack(I_output);')

catch

    if ~isempty(findstr(lasterr,' MEMORY '))
        disp('Insufficient Memory; See MATLAB HELP or try splitting input image into smaller segments')
    else
        disp(lasterr)
    end
end
