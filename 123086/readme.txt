This is the readme file for the model associated with the paper:

Carpenter, G.A. , Gaddam, C.S. , Mingolla, E., CONFIGR: A vision-
based model for long-range figure completion, Neural Networks -
Volume 20, 1109-1131 (2007).

Abstract

CONFIGR (CONtour FIgure GRound) is a computational model based on
principles of biological vision that completes sparse and noisy image
figures. Within an integrated vision/recognition system, CONFIGR
posits an initial recognition stage which identifies figure pixels
from spatially local input information. The resulting, and typically
incomplete, figure is fed back to the ?early vision? stage for long-
range completion via filling-in. The reconstructed image is then re-
presented to the recognition system for global functions such as
object recognition. In the CONFIGR algorithm, the smallest
independent image unit is the visible pixel, whose size defines a
computational spatial scale. Once the pixel size is fixed, the entire
algorithm is fully determined, with no additional parameter choices.
Multi-scale simulations illustrate the vision/recognition system.
Open-source CONFIGR code is available online, but all examples can be
derived analytically, and the design principles applied at each step
are transparent. The model balances filling-in as figure against
complementary filling-in as ground, which blocks spurious figure
completions. Lobe computations occur on a subpixel spatial scale.
Originally designed to fill-in missing contours in an incomplete
image such as a dashed line, the same CONFIGR system connects and
segments sparse dots, and unifies occluded objects from pieces
locally identified as figure in the initial recognition stage. The
model self-scales its completion distances, filling-in across gaps of
any length, where unimpeded, while limiting connections among dense
image-figure pixel groups that already have intrinsic form. Long-
range image completion promises to play an important role in adaptive
processors that reconstruct images from highly compressed video and
still camera images.

This software has been realized by Sai Gaddam in collaboration with
Gail Carpenter at the CNS Tehcnology Lab at Boston University
(http://techlab.bu.edu/).  Please see the license.txt file.

Software Description:

CONFIGR (CONtour FIgure and GRound) is a model that performs long-
range contour completion on large-scale images. CONFIGR accomplishes
this through a mechanism that fills-in both figure and ground via
complementary process.

Code Description

MATLAB 7.0 or higher. Usage: I_output = runCONFIGR(I,PixRes) This
runs CONFIGR with the following defaults: PixRes: pixel resolution
default=1: CONFIGR pixel resolution is the same as that of the input
image Advanced Options: 
I_output = runCONFIGR(I,PixRes,NumIter,ShrinkFact)
I: input image PixRes: pixel
resolution 1: fine 2: medium 3: coarse The following options provide
computational flexibility but are not model parameters. NumIter:
number of iterations CONFIGR simulation can be forced to stop early
by setting a low number of iterations. ShrinkFact: ratio of desired
image size to actual size Sparse images can be resized for faster
runtimes. Raw CONFIGR output (ground and figure-filled rectangles,
and interpolating diagonals) can be obtained using [I_output,
I_output_raw, Idiagonals]=runCONFIGR(I)

*CONFIGR (COntour FIgure and GRound)

*The following files are required to run CONFIGR:


CheckBound
CONFIGR_6.m
CONFIGR_6_CreatePlots.m
CONFIGR_6_FindBound.m
CONFIGR_6_Init.m
CONFIGR_6_LobePropagate.m
CreateNewDiagonals.m
EmptyRectangleTypeOne_Ground_6.m
EmptyRectangleTypeTwo_Ground_6.m
FillingFIGURE.m
FillingGROUND.m
plotIpixels.m
runCONFIGR.m
sign01.m
sign02.m
SimpleRoadImage.m

*Images are stored in the following .mat files

CONFIGRtestimages.mat
MontereyImage.mat
 
The code was tested using Matlab V 6.5 and v 7.5

*Comments
Images in the test folder provided here can be used as input
according to the following notation:



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
% I_output = runCONFIGR(I,PixRes)
% I ---------- input image
% PixRes ----- pixel resolution 1:fine 2:medium 3:coarse
%              N implies that each coarse-scale pixel is created
%              from an
%              2^Nx2^N square of fine-scale (original) pixels
% This runs CONFIGR with the following defaults:
% PixRes: pixel resolution default=1: CONFIGR pixel resolution is
% the same as that of the input image
%
%%%%%%%OPTIONS:%%%%%%%%%%%%%%%%
% I_output = runCONFIGR(I,PixRes,NumIter,ShrinkFact)
% I ----------- input image
% PixRes ------ pixel resolution 1:fine 2:medium 3:coarse
% NumIter ----- number of iterations CONFIGR simulation can be
% forced to stop early by setting a low number of iterations.
% ShrinkFact -- ratio of desired image size to actual size Sparse
% images can be resized for faster runtimes.
% Raw CONFIGR output (with ground and figure-filled rectangles)
% can be obtained using 
% [I_output, I_output_raw,Idiagonals]=runCONFIGR(I)
