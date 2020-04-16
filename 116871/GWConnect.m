function GWConnect
% A global workspace with an internal loop using spiking neurons
%
% Murray Shanahan
% July 2006 - September 2006
% With corrections January 2009
%
% GWConnect sets up the initial connections between all the areas in the
% model. Paths between areas are predefined. Paths within cortical areas
% are determined by STDP. The resulting network is saved in a file.


close all;

N = 16;     % Workspace areas are N by N matrices
Nc = 32;    % Cortical columns are Nc by Nc matrices
Dmax = 40;  % Maximum propagation delay

% Read in training stimuli
Image1 = imread('Blob1NW.bmp');
Image2 = imread('Blob2NE.bmp');
Image3 = imread('Blob1NE.bmp');
Image4 = imread('Blob2SE.bmp');
Image5 = imread('Blob2SW.bmp');

% Normalise pixels (stimuli are 16 colour bitmaps)
Image1 = double(Image1) ./ 15;
Image2 = double(Image2) ./ 15;
Image3 = double(Image3) ./ 15;
Image4 = double(Image4) ./ 15;
Image5 = double(Image5) ./ 15;

% Invert y axes for proper plotting
Image1(1:Nc,:) = Image1(Nc:-1:1,:);
Image2(1:Nc,:) = Image2(Nc:-1:1,:);
Image3(1:Nc,:) = Image3(Nc:-1:1,:);
Image4(1:Nc,:) = Image4(Nc:-1:1,:);
Image5(1:Nc,:) = Image5(Nc:-1:1,:);


% Neuron parameters: workspace areas

% Area W1
r = rand(N,N);
layer{1}.name = 'W1 (excitatory)';
layer{1}.pos = 6;
layer{1}.rows = N;
layer{1}.columns = N;
layer{1}.a = 0.02*ones(N,N);
layer{1}.b = 0.2*ones(N,N);
layer{1}.c = -65+15*r.^2;
layer{1}.d = 8-6*r.^2;
% Area W2
r = rand(N,N);
layer{2}.name = 'W2 (excitatory)';
layer{2}.pos = 7;
layer{2}.rows = N;
layer{2}.columns = N;
layer{2}.a = 0.02*ones(N,N);
layer{2}.b = 0.2*ones(N,N);
layer{2}.c = -65+15*r.^2;
layer{2}.d = 8-6*r.^2;
% Area W3
r = rand(N,N);
layer{3}.name = 'W3 (excitatory)';
layer{3}.pos = 8;
layer{3}.rows = N;
layer{3}.columns = N;
layer{3}.a = 0.02*ones(N,N);
layer{3}.b = 0.2*ones(N,N);
layer{3}.c = -65+15*r.^2;
layer{3}.d = 8-6*r.^2;
% Area W4
r = rand(N,N);
layer{4}.name = 'W4 (excitatory)';
layer{4}.pos = 9;
layer{4}.rows = N;
layer{4}.columns = N;
layer{4}.a = 0.02*ones(N,N);
layer{4}.b = 0.2*ones(N,N);
layer{4}.c = -65+15*r.^2;
layer{4}.d = 8-6*r.^2;
% Area W5
r = rand(N,N);
layer{5}.name = 'W5 (excitatory)';
layer{5}.pos = 10;
layer{5}.rows = N;
layer{5}.columns = N;
layer{5}.a = 0.02*ones(N,N);
layer{5}.b = 0.2*ones(N,N);
layer{5}.c = -65+15*r.^2;
layer{5}.d = 8-6*r.^2;

% Neuron parameters: workspace inhibitory populations

Nwi = floor(N*N/5);
% Area I1
r = rand(Nwi,1);
layer{9}.name = 'W1 (inhibitory)';
layer{9}.pos = 1;
layer{9}.rows = Nwi;
layer{9}.columns = 1;
layer{9}.a = 0.02+0.08*r;
layer{9}.b = 0.25-0.05*r;
layer{9}.c = -65*ones(Nwi,1);
layer{9}.d = 2*ones(Nwi,1);
% Area I2
r = rand(Nwi,1);
layer{10}.name = 'W2 (inhibitory)';
layer{10}.pos = 2;
layer{10}.rows = Nwi;
layer{10}.columns = 1;
layer{10}.a = 0.02+0.08*r;
layer{10}.b = 0.25-0.05*r;
layer{10}.c = -65*ones(Nwi,1);
layer{10}.d = 2*ones(Nwi,1);
% Area I3
r = rand(Nwi,1);
layer{11}.name = 'W3 (inhibitory)';
layer{11}.pos = 3;
layer{11}.rows = Nwi;
layer{11}.columns = 1;
layer{11}.a = 0.02+0.08*r;
layer{11}.b = 0.25-0.05*r;
layer{11}.c = -65*ones(Nwi,1);
layer{11}.d = 2*ones(Nwi,1);
% Area I4
r = rand(Nwi,1);
layer{12}.name = 'W4 (inhibitory)';
layer{12}.pos = 4;
layer{12}.rows = Nwi;
layer{12}.columns = 1;
layer{12}.a = 0.02+0.08*r;
layer{12}.b = 0.25-0.05*r;
layer{12}.c = -65*ones(Nwi,1);
layer{12}.d = 2*ones(Nwi,1);
% Area I5
r = rand(Nwi,1);
layer{13}.name = 'W5 (inhibitory)';
layer{13}.pos = 5;
layer{13}.rows = Nwi;
layer{13}.columns = 1;
layer{13}.a = 0.02+0.08*r;
layer{13}.b = 0.25-0.05*r;
layer{13}.c = -65*ones(Nwi,1);
layer{13}.d = 2*ones(Nwi,1);

% Neuron parameters: cortical columns

Ni = floor(Nc/5);  % Ni by Nc Inhibitory neurons per column
Ne = Nc-Ni;        % Ne by Nc Excitatory neurons per column
% Column 1
re = rand(Ne,Nc);
ri = rand(Ni,Nc);
layer{6}.name = 'C1';
layer{6}.pos = 16;
layer{6}.rows = Nc;
layer{6}.columns = Nc;
layer{6}.a(1:Ne,1:Nc) = 0.02*ones(Ne,Nc);
layer{6}.b(1:Ne,1:Nc) = 0.2*ones(Ne,Nc);
layer{6}.c(1:Ne,1:Nc) = -65+15*re.^2;
layer{6}.d(1:Ne,1:Nc) = 8-6*re.^2;
layer{6}.a(Ne+1:Nc,1:Nc) = 0.02+0.08*ri;
layer{6}.b(Ne+1:Nc,1:Nc) = 0.25-0.05*ri;
layer{6}.c(Ne+1:Nc,1:Nc) = -65*ones(Ni,Nc);
layer{6}.d(Ne+1:Nc,1:Nc) = 2*ones(Ni,Nc);
% Column 2
re = rand(Ne,Nc);
ri = rand(Ni,Nc);
layer{7}.name = 'C2';
layer{7}.pos = 17;
layer{7}.rows = Nc;
layer{7}.columns = Nc;
layer{7}.a(1:Ne,1:Nc) = 0.02*ones(Ne,Nc);
layer{7}.b(1:Ne,1:Nc) = 0.2*ones(Ne,Nc);
layer{7}.c(1:Ne,1:Nc) = -65+15*re.^2;
layer{7}.d(1:Ne,1:Nc) = 8-6*re.^2;
layer{7}.a(Ne+1:Nc,1:Nc) = 0.02+0.08*ri;
layer{7}.b(Ne+1:Nc,1:Nc) = 0.25-0.05*ri;
layer{7}.c(Ne+1:Nc,1:Nc) = -65*ones(Ni,Nc);
layer{7}.d(Ne+1:Nc,1:Nc) = 2*ones(Ni,Nc);
% Column 3
re = rand(Ne,Nc);
ri = rand(Ni,Nc);
layer{8}.name = 'C3';
layer{8}.pos = 18;
layer{8}.rows = Nc;
layer{8}.columns = Nc;
layer{8}.a(1:Ne,1:Nc) = 0.02*ones(Ne,Nc);
layer{8}.b(1:Ne,1:Nc) = 0.2*ones(Ne,Nc);
layer{8}.c(1:Ne,1:Nc) = -65+15*re.^2;
layer{8}.d(1:Ne,1:Nc) = 8-6*re.^2;
layer{8}.a(Ne+1:Nc,1:Nc) = 0.02+0.08*ri;
layer{8}.b(Ne+1:Nc,1:Nc) = 0.25-0.05*ri;
layer{8}.c(Ne+1:Nc,1:Nc) = -65*ones(Ni,Nc);
layer{8}.d(Ne+1:Nc,1:Nc) = 2*ones(Ni,Nc);

% Neuron parameters: lateral inhibition

Nci = floor(Nc*Nc/5);

% Lateral inhibition for column 2
r = rand(Nci,1);
layer{17}.name = 'L2';
layer{17}.pos = 22;
layer{17}.rows = Nci;
layer{17}.columns = 1;
layer{17}.a = 0.02+0.08*r;
layer{17}.b = 0.25-0.05*r;
layer{17}.c = -65*ones(Nci,1);
layer{17}.d = 2*ones(Nci,1);
% Lateral inhibition for column 3
r = rand(Nci,1);
layer{18}.name = 'L3';
layer{18}.pos = 23;
layer{18}.rows = Nci;
layer{18}.columns = 1;
layer{18}.a = 0.02+0.08*r;
layer{18}.b = 0.25-0.05*r;
layer{18}.c = -65*ones(Nci,1);
layer{18}.d = 2*ones(Nci,1);

% Neuron parameters: workspace access areas

% Access area 1
r = rand(N,N);
layer{14}.name = 'A1';
layer{14}.pos = 11;
layer{14}.rows = N;
layer{14}.columns = N;
layer{14}.a = 0.02*ones(N,N);
layer{14}.b = 0.2*ones(N,N);
layer{14}.c = -65+15*r.^2;
layer{14}.d = 8-6*r.^2;
% Access area 2
r = rand(N,N);
layer{15}.name = 'A2';
layer{15}.pos = 12;
layer{15}.rows = N;
layer{15}.columns = N;
layer{15}.a = 0.02*ones(N,N);
layer{15}.b = 0.2*ones(N,N);
layer{15}.c = -65+15*r.^2;
layer{15}.d = 8-6*r.^2;
% Access area 3
r = rand(N,N);
layer{16}.name = 'A3';
layer{16}.pos = 13;
layer{16}.rows = N;
layer{16}.columns = N;
layer{16}.a = 0.02*ones(N,N);
layer{16}.b = 0.2*ones(N,N);
layer{16}.c = -65+15*r.^2;
layer{16}.d = 8-6*r.^2;


% Connectivity matrices (synaptic weights)

% layer{i}.S{j} is the connectivity matrix from layer j to layer i
% S(i,j) is the strength of the connection from neuron j to neuron i

% Clear connectivity matrices
L = length(layer);
for i=1:L
   for j=1:L
      layer{i}.S{j} = [];
      layer{i}.factor{j} = [];
      layer{i}.delay{j} = [];
   end
end

% Paths between workspace areas
layer{4}.S{1} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{1}.S{5} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{5}.S{2} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{2}.S{1} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{1}.S{3} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{3}.S{2} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{2}.S{4} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{4}.S{3} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{3}.S{5} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{5}.S{4} = spdiags(ones(N*N,1),0,N*N,N*N);

% Workspace inhibitory paths
layer{9}.S{1} = 1.0*rand(Nwi,N*N);
layer{4}.S{9} = -1.0*rand(N*N,Nwi);
layer{2}.S{9} = -1.0*rand(N*N,Nwi);
layer{10}.S{2} = 1.0*rand(Nwi,N*N);
layer{5}.S{10} = -1.0*rand(N*N,Nwi);
layer{3}.S{10} = -1.0*rand(N*N,Nwi);
layer{11}.S{3} = 1.0*rand(Nwi,N*N);
layer{1}.S{11} = -1.0*rand(N*N,Nwi);
layer{4}.S{11} = -1.0*rand(N*N,Nwi);
layer{12}.S{4} = 1.0*rand(Nwi,N*N);
layer{2}.S{12} = -1.0*rand(N*N,Nwi);
layer{5}.S{12} = -1.0*rand(N*N,Nwi);
layer{13}.S{5} = 1.0*rand(Nwi,N*N);
layer{3}.S{13} = -1.0*rand(N*N,Nwi);
layer{1}.S{13} = -1.0*rand(N*N,Nwi);

% Recurrent connections in cortical columns
we = 0.1;
wi = -0.5;
M = Nc-floor(Nc/5);
s1 = zeros(Nc*Nc,Nc*Nc);
s2 = zeros(Nc*Nc,Nc*Nc);
s3 = zeros(Nc*Nc,Nc*Nc);
% 80% excitatory neurons - rows 1 to M of an Nc by Nc map
for i=1:M
   for j=1:Nc
      s1(:,Nc*(j-1)+i) = we*rand(Nc*Nc,1);
      s2(:,Nc*(j-1)+i) = we*rand(Nc*Nc,1);
      s3(:,Nc*(j-1)+i) = we*rand(Nc*Nc,1);
   end
end
% 20% inhibitory neurons - rows M+1 to Nc of an Nc by Nc map
for i=M+1:Nc
   for j=1:Nc
      s1(:,Nc*(j-1)+i) = wi*rand(Nc*Nc,1);
      s2(:,Nc*(j-1)+i) = wi*rand(Nc*Nc,1);
      s3(:,Nc*(j-1)+i) = wi*rand(Nc*Nc,1);
   end
end
layer{6}.S{6} = s1;
layer{7}.S{7} = s2;
layer{8}.S{8} = s3;

% Lateral inhibition
s1 = sparse(Nci,Nc*Nc);
s2 = sparse(Nci,Nc*Nc);
for i=1:N
   for j=1:N
      s1(:,Nc*(i-1+Nc/2)+j) = 1.0*rand(Nci,1);
      s2(:,Nc*(i-1+Nc/2)+j) = 1.0*rand(Nci,1);
   end
end
layer{17}.S{7} = s1;
layer{8}.S{17} = -1.0*rand(Nc*Nc,Nci);
layer{18}.S{8} = s2;
layer{7}.S{18} = -1.0*rand(Nc*Nc,Nci);
layer{17}.S{18} = -1.0*rand(Nci,Nci);
layer{18}.S{17} = -1.0*rand(Nci,Nci);

% Top-down amplification
s1 = sparse(Nc*Nc,N*N);
for i=1:N
   for j=1:N
      s1(Nc*(i-1+Nc/2)+j,N*(i-1)+j) = 1;
   end
end
layer{6}.S{14} = s1;
layer{7}.S{15} = s1;
layer{8}.S{16} = s1;

% Paths to and from cortical columns
s1 = sparse(Nc*Nc,N*N);
s2 = sparse(N*N,Nc*Nc);
for i=1:N
   for j=1:N
      s1(Nc*(i-1)+j,N*(i-1)+j) = 1;
      s2(N*(i-1)+j,Nc*(i-1+Nc/2)+j) = 1;
   end
end
layer{6}.S{1} = s1;
layer{14}.S{6} = s2;
layer{7}.S{2} = s1;
layer{15}.S{7} = s2;
layer{8}.S{2} = s1;
layer{16}.S{8} = s2;

% Paths from access areas into workspace
layer{1}.S{14} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{2}.S{15} = spdiags(ones(N*N,1),0,N*N,N*N);
layer{2}.S{16} = spdiags(ones(N*N,1),0,N*N,N*N);


% Connectivity scaling factors - Supplementary Note values

Fwc = 40;
Faw = 90;
Fww = 90;
Fcc = 3.0;
Fcl = 30;
Flc = 300;
Fwi = 4.5;
Fiw = 4.5;
Fac = 90;
Fll = 90;
Fca = 90;

% Paths between workspace areas
layer{4}.factor{1} = Fww;
layer{1}.factor{5} = Fww;
layer{5}.factor{2} = Fww;
layer{2}.factor{1} = Fww;
layer{1}.factor{3} = Fww;
layer{3}.factor{2} = Fww;
layer{2}.factor{4} = Fww;
layer{4}.factor{3} = Fww;
layer{3}.factor{5} = Fww;
layer{5}.factor{4} = Fww;

% Workspace inhibitory paths
layer{9}.factor{1} = Fwi;
layer{4}.factor{9} = Fiw;
layer{2}.factor{9} = Fiw;
layer{10}.factor{2} = Fwi;
layer{5}.factor{10} = Fiw;
layer{3}.factor{10} = Fiw;
layer{11}.factor{3} = Fwi;
layer{1}.factor{11} = Fiw;
layer{4}.factor{11} = Fiw;
layer{12}.factor{4} = Fwi;
layer{2}.factor{12} = Fiw;
layer{5}.factor{12} = Fiw;
layer{13}.factor{5} = Fwi;
layer{3}.factor{13} = Fiw;
layer{1}.factor{13} = Fiw;

% Paths to, from, and within cortical columns
layer{6}.factor{1} = Fwc;
layer{14}.factor{6} = Fca;
layer{6}.factor{6} = Fcc;
layer{7}.factor{2} = Fwc;
layer{15}.factor{7} = Fca;
layer{7}.factor{7} = Fcc;
layer{8}.factor{2} = Fwc;
layer{16}.factor{8} = Fca;
layer{8}.factor{8} = Fcc;

% Lateral inhibition paths
layer{17}.factor{7} = Fcl;
layer{8}.factor{17} = Flc;
layer{18}.factor{8} = Fcl;
layer{7}.factor{18} = Flc;
layer{17}.factor{18} = Fll;
layer{18}.factor{17} = Fll;

% Top-down amplification
layer{6}.factor{14} = Fac;
layer{7}.factor{15} = Fac;
layer{8}.factor{16} = Fac;

% Paths from access areas into workspace
layer{1}.factor{14} = Faw;
layer{2}.factor{15} = Faw;
layer{2}.factor{16} = Faw;


% Propagation delays

Daw = 20;
Dwc = 10;
Dww1 = 5;
Dww2 = 6;
Diw1 = 5;
Diw2 = 6;
Dac = 2;
Dcl = 2;
Dlc = 2;
Dll = 2;
Dwi = 2;
Dca = 2;

% Paths between workspace areas
layer{4}.delay{1} = ones(N*N,N*N)*Dww2;
layer{1}.delay{5} = ones(N*N,N*N)*Dww1;
layer{5}.delay{2} = ones(N*N,N*N)*Dww2;
layer{2}.delay{1} = ones(N*N,N*N)*Dww1;
layer{1}.delay{3} = ones(N*N,N*N)*Dww2;
layer{3}.delay{2} = ones(N*N,N*N)*Dww1;
layer{2}.delay{4} = ones(N*N,N*N)*Dww2;
layer{4}.delay{3} = ones(N*N,N*N)*Dww1;
layer{3}.delay{5} = ones(N*N,N*N)*Dww2;
layer{5}.delay{4} = ones(N*N,N*N)*Dww1;

% Workspace inhibitory populations
layer{9}.delay{1} = ones(Nwi,N*N)*Dwi;
layer{4}.delay{9} = ones(N*N,Nwi)*Diw2;
layer{2}.delay{9} = ones(N*N,Nwi)*Diw1;
layer{10}.delay{2} = ones(Nwi,N*N)*Dwi;
layer{5}.delay{10} = ones(N*N,Nwi)*Diw2;
layer{3}.delay{10} = ones(N*N,Nwi)*Diw1;
layer{11}.delay{3} = ones(Nwi,N*N)*Dwi;
layer{1}.delay{11} = ones(N*N,Nwi)*Diw2;
layer{4}.delay{11} = ones(N*N,Nwi)*Diw1;
layer{12}.delay{4} = ones(Nwi,N*N)*Dwi;
layer{2}.delay{12} = ones(N*N,Nwi)*Diw2;
layer{5}.delay{12} = ones(N*N,Nwi)*Diw1;
layer{13}.delay{5} = ones(Nwi,N*N)*Dwi;
layer{3}.delay{13} = ones(N*N,Nwi)*Diw2;
layer{1}.delay{13} = ones(N*N,Nwi)*Diw1;

% Paths to, from, and within cortical columns
layer{6}.delay{1} = ones(Nc*Nc,N*N)*Dwc;
layer{14}.delay{6} = ones(N*N,Nc*Nc)*Dca;
layer{6}.delay{6} = ceil(rand(Nc*Nc,Nc*Nc)*Dmax);
layer{7}.delay{2} = ones(Nc*Nc,N*N)*Dwc;
layer{15}.delay{7} = ones(N*N,Nc*Nc)*Dca;
layer{7}.delay{7} = ceil(rand(Nc*Nc,Nc*Nc)*Dmax);
layer{8}.delay{2} = ones(Nc*Nc,N*N)*Dwc;
layer{16}.delay{8} = ones(N*N,Nc*Nc)*Dca;
layer{8}.delay{8} = ceil(rand(Nc*Nc,Nc*Nc)*Dmax);

% Lateral inhibition paths
layer{17}.delay{7} = ones(Nci,Nc*Nc)*Dcl;
layer{8}.delay{17} = ones(Nc*Nc,Nci)*Dlc;
layer{18}.delay{8} = ones(Nci,Nc*Nc)*Dcl;
layer{7}.delay{18} = ones(Nc*Nc,Nci)*Dlc;
layer{17}.delay{18} = ones(Nci,Nci)*Dll;
layer{18}.delay{17} = ones(Nci,Nci)*Dll;

% Top-down amplification
layer{6}.delay{14} = ones(Nc*Nc,N*N)*Dac;
layer{7}.delay{15} = ones(Nc*Nc,N*N)*Dac;
layer{8}.delay{16} = ones(Nc*Nc,N*N)*Dac;

% Paths from access areas into workspace
layer{1}.delay{14} = ones(N*N,N*N)*Daw;
layer{2}.delay{15} = ones(N*N,N*N)*Daw;
layer{2}.delay{16} = ones(N*N,N*N)*Daw;


Icpulse = 10; % pulse current


% TRAIN

Fin = 40; % scaling factor for thalamic input

for r=1:2
   
   % Initial membrane potentials
   for lr=1:length(layer)
      layer{lr}.v = -65*ones(layer{lr}.rows,layer{lr}.columns);
      layer{lr}.u = layer{lr}.b.*layer{lr}.v;
   end

   % Clear lists of spikes
   for lr=1:length(layer)
      layer{lr}.firings = [];
   end
   
   for t=1:200

      % Display time every 10ms
      if mod(t,10) == 0
         t
      end

      % Column 1
      layer{6}.I = 1*randn(Nc,Nc);
      if (t == 20 || t == 25 || t == 30 || t == 35)
         layer{6}.I = layer{6}.I+Image1*Icpulse*Fin;
      end
      if (t == 80 || t == 85 || t == 90 || t == 95)
         layer{6}.I = layer{6}.I+Image2*Icpulse*Fin;
      end

      % Column 2
      layer{7}.I = 1*randn(Nc,Nc);
      if (t == 20 || t == 25 || t == 30 || t == 35)
         layer{7}.I = layer{7}.I+Image3*Icpulse*Fin;
      end
      if (t == 80 || t == 85 || t == 90 || t == 95)
         layer{7}.I = layer{7}.I+Image4*Icpulse*Fin;
      end
   
      % Column 3
      layer{8}.I = 1*randn(Nc,Nc);
      if (t == 20 || t == 25 || t == 30 || t == 35)
         layer{8}.I = layer{8}.I+Image3*Icpulse*Fin;
      end
      if (t == 80 || t == 85 || t == 90 || t == 95)
         layer{8}.I = layer{8}.I+Image5*Icpulse*Fin;
      end

      for lr=[6 7 8];
         layer = update(layer,lr,t,Dmax);
         % plot_layer(layer,lr,t);
         layer = STDPtrain(layer,lr,t);
      end

      % waitforbuttonpress;
      drawnow;
      
   end
   
end

save('Network.mat','layer');

end