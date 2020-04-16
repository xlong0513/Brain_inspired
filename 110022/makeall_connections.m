% makeall_connections adds all connections to the FEF local circuit.
% All parameters are defined within this file.
% 
% created: Jakob Heinzle 01/07

% General parameters
w_dist=1;  % distribution of weights w_dist=0:no distribution, 1:uniform, 2:gaussian (not yet implemented)
sigma_dist=0.5; % sigma of distribution, needs to be <=1
disp('Weight matrices are assigned.');

%----------------------------------------------------------------------
% to layer 4
%----------------------------------------------------------------------
% from E4 to E4
type='exc';from='E4';to='E4';weight=0.008;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0)+0.05*(diag(ones(nfac-1,1),1)+diag(ones(nfac-1,1),-1));
add_connection;
% from I4 to E4
type='inh';from='I4';to='E4';weight=0.06;tau=3;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from E6A to E4
type='exc';from='E6A';to='E4';weight=0.004;tau=5;sparseness=0.5;
small_matrix=antidiagonal_block(nfac,nfac,nfac);
add_connection;

% from E4 to I4
type='exc';from='E4';to='I4';weight=0.0025;tau=5;sparseness=0.25;
small_matrix=ones(nfac,nfac);
add_connection;
type='exc';from='E6S';to='I4';weight=0.004;tau=10;sparseness=0.5;
small_matrix=ones(nfac,nfac);
add_connection;
type='exc';from='E6S';to='I4';weight=0.0008;tau=50;sparseness=0.5;
small_matrix=antidiagonal_block(nfac,nfac,nfac);
add_connection;
type='exc';from='E23';to='I4';weight=0.001;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0)+diag(ones(nfac-1,1),1)+diag(ones(nfac-1,1),-1);
add_connection;

%----------------------------------------------------------------------
% to layer 2/3
%----------------------------------------------------------------------
% from E23 to E23
type='exc';from='E23';to='E23';weight=0.0048;tau=10;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from I4 to E4
type='inh';from='I23';to='E23';weight=0.08;tau=3;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from E4 to E23
type='exc';from='E4';to='E23';weight=0.0016;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from E5B to E23
type='exc';from='E5B';to='E23';weight=0.0085;tau=10;sparseness=0.5;
small_matrix=zeros(nfac,nfac);small_matrix((nfac+1)/2,:)=1;
add_connection;
% from E4 to I4
type='exc';from='E23';to='I23';weight=0.002;tau=5;sparseness=0.25;
small_matrix=ones(nfac,nfac);
add_connection;
% from E5B to I23
type='exc';from='E5B';to='I23';weight=0.02;tau=5;sparseness=0.5;
small_matrix=ones(nfac,nfac);small_matrix((nfac+1)/2,:)=0;
add_connection;
% from ERb to I23
type='exc';from='ERb';to='I23';weight=0.008;tau=10;sparseness=0.5;
small_matrix=diag(ones(nfac,1));
add_connection;

%----------------------------------------------------------------------
% to layer 5
%--------------------------------------------------------------------------
% from E5R to E5R
type='exc';from='E5R';to='E5R';weight=0.002;tau=50;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
type='exc';from='E23';to='E5R';weight=0.00132;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);small_matrix((nfac+1)/2,(nfac+1)/2)=0;
add_connection;
% from IFIX to E5R
type='inh';from='IFIX';to='E5R';weight=0.0035;tau=3;sparseness=0.5;
small_matrix=ones(nfac,1);
add_connection;
% from I5B to E5R
type='inh';from='I5B';to='E5R';weight=0.02;tau=10;sparseness=0.5;
small_matrix=diag(ones(nfac,1));
add_connection;
% from E5R to I5R
type='exc';from='E5R';to='I5R';weight=0.015;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1));
add_connection;

% from E5B to E5B
type='exc';from='E5B';to='E5B';weight=0.06;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from E5R to E5B
type='exc';from='E5R';to='E5B';weight=0.01;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from I5B to E5B
type='inh';from='I5B';to='E5B';weight=0.125;tau=3;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from E5B to I5B
type='exc';from='E5B';to='I5B';weight=0.05;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;

%----------------------------------------------------------------------
% to layer 6
%-------------------------------------------------------------------------
% from E5B to E6S
type='exc';from='E5B';to='E6S';weight=0.04;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1));
add_connection;
% from E23 to E6A
type='exc';from='E23';to='E6A';weight=0.005;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1));
add_connection;
%from EFa to E6A
type='exc';from='EFa';to='E6A';weight=0.008;tau=5;sparseness=0.5;
small_matrix=ones(nfac,nfac);
add_connection;
%from EFf to E6A
type='exc';from='EFf';to='E6A';weight=0.014;tau=5;sparseness=0.5;
small_matrix=zeros(nfac,nfac);small_matrix((nfac+1)/2,:)=1;
add_connection;

%----------------------------------------------------------------------
% to fixation neurons
%--------------------------------------------------------------------------
% from I5R to IFIX
type='inh';from='I5R';to='IFIX';weight=0.05;tau=3;sparseness=0.5;
small_matrix=ones(1,nfac);
add_connection;
% from E23 to IFIX
type='exc';from='E23';to='IFIX';weight=0.002;tau=5;sparseness=0.5;
small_matrix=zeros(1,nfac);small_matrix((nfac+1)/2)=1;
add_connection;

%----------------------------------------------------------------------
% Connections of the Recognition Module
%----------------------------------------------------------------------

% from EFp to ERr
type='exc';from='EFp';to='ERr';weight=0.002;tau=5;sparseness=0.5;
small_matrix=zeros(nfac,nfac);
small_matrix((nfac+1)/2,(nfac+1)/2)=1;
add_connection;
% from EFf to ERr
type='exc';from='EFf';to='ERr';weight=0.002;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1));
add_connection;
% from EFa to ERr
type='exc';from='EFa';to='ERr';weight=0.002;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1));
add_connection;
% from EFspace to ERr
type='exc';from='EFspace';to='ERr';weight=0.0025;tau=5;sparseness=0.5;
small_matrix=zeros(nfac,1);small_matrix((nfac+1)/2)=1;
add_connection;

%from IF to EFp
type='inh';from='IF';to='EFp';weight=0.3;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
%from IF to EFf
type='inh';from='IF';to='EFf';weight=0.3;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
%from IF to EFa
type='inh';from='IF';to='EFa';weight=0.3;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
%from IF to EFspace
type='inh';from='IF';to='EFspace';weight=0.3;tau=5;sparseness=0.5;
small_matrix=zeros(1,nfac); small_matrix((nfac+1)/2)=1;
add_connection;

% from E23 to IF
type='inh';from='E23';to='IF';weight=0.012;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;

% from ERr to ERr
type='exc';from='ERr';to='ERr';weight=0.0006;tau=50;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from IRb to ERr
type='inh';from='IRb';to='ERr';weight=0.12;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;


% from ERb to ERb
type='exc';from='ERb';to='ERb';weight=0.014;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from IRb to ERb
type='inh';from='IRb';to='ERb';weight=0.04;tau=3;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from ERr to ERb
type='exc';from='ERr';to='ERb';weight=0.006;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;
% from ERb to IRb
type='exc';from='ERb';to='IRb';weight=0.02;tau=5;sparseness=0.5;
small_matrix=diag(ones(nfac,1),0);
add_connection;