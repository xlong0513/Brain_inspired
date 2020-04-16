% spike_to_rate takes the data from simulations where the spiking data was
% saved  and transforms the spiking data to populations rates.
% These rates are stored as matrices called ..._RATE. The size of the matrix 
% is (npop x simulationtime). The plot command plot(..._RATE') will plot all 
% population rates of a single layer over the whole time period.
%
% created: Jakob Heinzle 01/07

clear *RATE

nfac=21;
foveaplot=1;
T=length(ERr(1,:));
poolsize=100;
factE=1000/poolsize;
factI=1000/poolsize*4;
factE5=factE*5/2;
factIFIX=factE;

% compute the population firing of all the populations.
E4f=full(sparse(ceil(E4(:,2)/poolsize),ceil(E4(:,1)-0.01),1,nfac,T));
E23f=full(sparse(ceil(E23(:,2)/poolsize),ceil(E23(:,1)-0.01),1,nfac,T));
E5Rf=full(sparse(ceil(E5R(:,2)/poolsize*5/2),ceil(E5R(:,1)-0.01),1,nfac,T));
E5Bf=full(sparse(ceil(E5B(:,2)/poolsize*5/2),ceil(E5B(:,1)-0.01),1,nfac,T));
E6Af=full(sparse(ceil(E6A(:,2)/poolsize*2),ceil(E6A(:,1)-0.01),1,nfac,T));
E6Sf=full(sparse(ceil(E6S(:,2)/poolsize/2),ceil(E6S(:,1)-0.01),1,nfac,T));

I4f=full(sparse(ceil(I4(:,2)/poolsize*4),ceil(I4(:,1)-0.01),1,nfac,T));
I23f=full(sparse(ceil(I23(:,2)/poolsize*4),ceil(I23(:,1)-0.01),1,nfac,T));
I5Rf=full(sparse(ceil(I5R(:,2)/poolsize*4),ceil(I5R(:,1)-0.01),1,nfac,T));
I5Bf=full(sparse(ceil(I5B(:,2)/poolsize*4),ceil(I5B(:,1)-0.01),1,nfac,T));

IFIXf=full(sparse(ceil(IFIX(:,2)/poolsize),ceil(IFIX(:,1)-0.01),1,1,T));

for i=1:nfac
    I23_RATE(i,:)=smooth_synapse(I23f(i,:),1,10)'*factI; E23_RATE(i,:)=smooth_synapse(E23f(i,:),1,10)'*factE;
    I4_RATE(i,:)=smooth_synapse(I4f(i,:),1,10)'*factI; E4_RATE(i,:)=smooth_synapse(E4f(i,:),1,10)'*factE;
    I5R_RATE(i,:)=smooth_synapse(I5Rf(i,:),1,10)'*factI; E5R_RATE(i,:)=smooth_synapse(E5Rf(i,:),1,10)'*factE5;
    I5B_RATE(i,:)=smooth_synapse(I5Bf(i,:),1,10)'*factI; E5B_RATE(i,:)=smooth_synapse(E5Bf(i,:),1,10)'*factE5;
    E6A_RATE(i,:)=smooth_synapse(E6Af(i,:),1,10)'*factE*2;E6S_RATE(i,:)=smooth_synapse(E6Sf(i,:),1,10)*factE*2;
    EFp_RATE(i,:)=smooth_synapse(EFp(i,:),1,10)'; EFf_RATE(i,:)=smooth_synapse(EFf(i,:),1,10)'; 
    EFa_RATE(i,:)=smooth_synapse(EFa(i,:),1,10)';
    IF_RATE(i,:)=smooth_synapse(IF(i,:),1,10)'; 
    ERr_RATE(i,:)=smooth_synapse(ERr(i,:),1,10)';ERb_RATE(i,:)=smooth_synapse(ERb(i,:),1,10)';
    IRb2_RATE(i,:)=smooth_synapse(IRb(i,:),1,10)';  
end
IFIX_RATE=smooth_synapse(IFIXf,1,10)'*factIFIX;
clear E4f I4f E23f I23f E5Rf I5Rf E5Bf I5Bf E6Af E6Sf IFIXf

