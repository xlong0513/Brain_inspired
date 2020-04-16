% ReadOutVariables creates the variables for the readout that are then
% saved.
%
% created: Jakob Heinzle 01/07

% create variables, where the spiking is saved.
E4=zeros(tmax*10,2);
I4=zeros(tmax*20,2);
E23=zeros(tmax*10,2);
I23=zeros(tmax*20,2);
E5B=zeros(tmax,2);
I5B=zeros(tmax,2);
E5R=zeros(tmax*10,2);
I5R=zeros(tmax*10,2);
E6A=zeros(tmax*10,2);
E6S=zeros(tmax*10,2);
IFIX=zeros(tmax*10,2);

% neurons in the REC module
EFp=zeros(nfac,tmax);
EFf=zeros(nfac,tmax);
EFa=zeros(nfac,tmax);
IF=zeros(nfac,tmax);
ERr=zeros(nfac,tmax);
ERb=zeros(nfac,tmax);
IRb=zeros(nfac,tmax);

% auxiliary rates used for plotting.
E4_HZ=zeros(nfac,tmax);
E23_HZ=zeros(nfac,tmax);
E5B_HZ=zeros(nfac,tmax);
E5R_HZ=zeros(nfac,tmax);