% Setup_Demo is used to set the task, the visual input and the plotting
% of the simulation. Only this file has to changed when running a
% simulation with FEF_DEMO.
%
% created: Jakob Heinzle 04/07

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select the task you want to run and define the input.
% 1=Pro-saccade, 2=Fixation, 3=Anti-Saccade, (reactive saccades)
% 4=Scanning, 
% 5=delayed memory pro-saccade, 6= delayed memory anti-saccade
% 7=delayed memory anti-saccade (task stimulus shown during delay)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

select_task=3; % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the input arrays and set the feature assigned to the target.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if select_task==1
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   disp('% Pro-saccade task');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   inarray1=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   featarray=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   tfix_on=0;
   tfix_off=50;
   tvisual_on=50;
   tvisual_off=250;
   tvisual_transoff=90;
   tmax=500;
elseif select_task==2
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   disp('% No-go task');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   inarray1=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   featarray=[0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   tfix_on=0;
   tfix_off=50;
   tvisual_on=50;
   tvisual_off=250;
   tvisual_transoff=90;
   tmax=500;
elseif select_task==3
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   disp('% Anti_Saccade task');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   inarray1=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   featarray=[0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   tfix_on=0;
   tfix_off=50;
   tvisual_on=50;
   tvisual_off=250;
   tvisual_transoff=90;
   tmax=500;
elseif select_task==4
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   disp('% Scanning task');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   inarray1=[0 0 0 0 0.9 0 1 0 0.8 0 0 0 1 0 0 0 0.8 0 0.9 0 0];
   featarray=(inarray1>0);
   tfix_on=0;
   tfix_off=50;
   tmax=10000;
   tvisual_on=50;
   tvisual_off=tmax;
   tvisual_transoff=90;        
elseif select_task==5
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   disp('% Delayed memory pro-saccade task');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   inarray1=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   featarray=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   tfix_on=0;
   tfix_off=650;
   tvisual_on=50;
   tvisual_off=250;
   tvisual_transoff=90;
   tmax=1000;
elseif select_task==6
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   disp('% Delayed memory anti-saccade task:');
   disp('% Task defined with first stimulus presentation');
   disp('% Compare paper by Zhang and Barash, 2004');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   inarray1=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   featarray=[0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   tfix_on=0;
   tfix_off=650;
   tvisual_on=50;
   tvisual_off=350;
   tvisual_transoff=90;
   tmax=1000;
elseif select_task==7
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   disp('% Delayed memory anti-saccade task');
   disp('% Task defined with second stimulus presentation');
   disp('% Compare paper by Amemori and Sawaguchi, 2006');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   inarray1=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   featarray=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
   tfix_on=0;
   tfix_off=950;
   tvisual_on=50;
   tvisual_off=250;
   tvisual_transoff=90;
   tmax=1200;
   stimdelay=400; % defines the delay (compared to the visual input,
   % with which the task stimulus is shown to the monkey.
end

%=====================================================================
% Printing options (can also be defined manually
% 
% The population given here will be plotted at the end of the simulations: 
% -10 is the most left position, 0 the fovea and 10 the most right
% position.
%=====================================================================
if select_task==4
    print_positions=-10:10; % plot all positions
else
   print_positions=[-4 0 4];  %only relevant positions for these tasks.
end

%=====================================================================
% Initialize graphics (displays only a single layer, here layer 4)
%=====================================================================

skip = 5;  % refresh interval for plot in ms;
t=0;
figure(1);
clf;
subplot(4,1,1);
hInp=bar(-10:10,zeros(1,21));
axis([-10.5 10.5 0 1.5]);
hold on;
xlabel('Retinotopic position')
ylabel('Input strength')
title('Input to layer 4')
hText=text(-10,1.1,[int2str(round(t)),' ms'],'FontSize',18);

subplot(4,1,2);
hL4=bar(-10:10,zeros(1,21));
axis([-10.5 10.5 0 100]);
hold on;
xlabel('Retinotopic position')
ylabel('Firing Rate')
title('Layer 4E')

subplot(4,1,3);
hL23=bar(-10:10,zeros(1,21));
axis([-10.5 10.5 0 100]);
hold on;
xlabel('Retinotopic position')
ylabel('Firing Rate')
title('Layer 2/3 E')

subplot(4,1,4);
hL5=bar(-10:10,zeros(1,21));
axis([-10.5 10.5 0 100]);
hold on;
xlabel('Retinotopic position')
ylabel('Firing Rate')
title('Layer 5 E')

%====================================================================
% Define which values to save after the simulation.
%====================================================================

savepath=['save ''FEF_RESULT',int2str(select_task),''' E4 I4 E23 I23 E5R E5B I5R I5B E6A E6S ERb ERr IRb IF EFp EFf EFf ATT IFIX FOVEA'];
