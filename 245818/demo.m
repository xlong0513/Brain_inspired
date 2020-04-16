%%% Example code for HV_disparity code
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


n_scales =9;  %pyramidal scales
th=1e-4;      %V1 global energy threshold 
th2=2.54;     %V1 energy threshold 
n_filters=18; %spatial filer orientations
Dori=2;       %MT (disparity) directions


load LeftRightImages00

D = HV_disparity(II,n_scales,n_filters,th,th2,Dori);

figure,imagesc(II(:,:,1)),colormap gray, title('left image')
figure,imagesc(D(:,:,1)),colormap jet, title('horizontal disparity (estimated)')
figure,imagesc(D(:,:,2)),colormap jet, title('vertical disparity (estimated)')

