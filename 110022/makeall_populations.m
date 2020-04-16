% makeall_populations defines all the background parameters and the size of 
% populations within for the simulation of the FEF network
%
% created: Jakob Heinzle 01/07

%-------------------------------------------------------------------------
% Define populations within the FEF.
%-------------------------------------------------------------------------
%background input feature selection neurons
gmaxE_ext = 0.02; %
gmaxI_ext = 3*gmaxE_ext;

% Add populations in Layer 4.
name='E4';type='exc';poolsize=100;nretpos=21;bgE=0.472;bgI=0.34;
add_population;
name='I4';type='inh';poolsize=25;nretpos=21;bgE=0.46;bgI=0.4;
add_population;

% Add populations in Layer 2/3.
name='E23';type='exc';poolsize=100;nretpos=21;bgE=0.472;bgI=0.34;
add_population;
name='I23';type='inh';poolsize=25;nretpos=21;bgE=0.46;bgI=0.4;
add_population;

% Add populations in Layer 5.
name='E5R';type='exc';poolsize=40;nretpos=21;bgE=0.45;bgI=0.34;
add_population;
name='I5R';type='inh';poolsize=25;nretpos=21;bgE=0.42;bgI=0.34;
add_population;

name='E5B';type='exc';poolsize=40;nretpos=21;bgE=0.38;bgI=0.30;
add_population;
name='I5B';type='inh';poolsize=25;nretpos=21;bgE=0.32;bgI=0.34;
add_population;

% Add populations in Layer 6.
name='E6S';type='exc';poolsize=50;nretpos=21;bgE=0.44;bgI=0.34;
add_population;
name='E6A';type='exc';poolsize=50;nretpos=21;bgE=0.2;bgI=0.34;
add_population;

% Add Fixation neurons
name='IFIX';type='inh';poolsize=100;nretpos=1;bgE=0.46;bgI=0.12;
add_population;


%-------------------------------------------------------------------------
% define populations of the Recognition Module.
%-------------------------------------------------------------------------

% Add inhibory neurons in feature detection
name='IF';type='inh';poolsize=25;nretpos=21;bgE=0.55;bgI=0.34;
add_population;

% Add feature detection arrays.
name='EFp';type='inh';poolsize=100;nretpos=21;bgE=0.42;bgI=0.30; %Pro-Saccade
add_population;
name='EFf';type='inh';poolsize=100;nretpos=21;bgE=0.42;bgI=0.30; %Fixation
add_population;
name='EFa';type='inh';poolsize=100;nretpos=21;bgE=0.42;bgI=0.30; %Anti-Saccade
add_population;
name='EFspace';type='inh';poolsize=100;nretpos=1;bgE=0.42;bgI=0.30; %Recognizing spaces
add_population;

% Add recognition arrays
name='ERr';type='exc';poolsize=100;nretpos=21;bgE=0.45;bgI=0.33;
add_population;
name='ERb';type='exc';poolsize=100;nretpos=21;bgE=0.38;bgI=0.30;
add_population;
name='IRb';type='inh';poolsize=25;nretpos=21;bgE=0.32;bgI=0.34;
add_population;
