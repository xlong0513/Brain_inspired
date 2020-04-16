% makeall_inputs defines all the inputs for the simulation of the FEF.
% These inputs are: Dorsal visual input and Fixation input to FEF; 
% Ventral feature specific input to Recognition Module.
%
% created: Jakob Heinzle 01/07


clear inputs
inputs.ninputs=0;

% Add the visual input to layer 4
name='VisualInput';type='exc';to='E4';inarray=inarray1;retinotopic=1;
MeanInp=0.056;NoiseLevel=1;t_on=tvisual_on;t_trans_off=tvisual_transoff;t_off=tvisual_off;sustained_level=0.5;
add_input;

% Add the fixation input to fixation neurons
name='FixationInput';type='exc';to='IFIX';inarray=1;retinotopic=0;
MeanInp=0.2;NoiseLevel=1;t_on=tfix_on;t_trans_off=tmax;t_off=tfix_off;sustained_level=1;
add_input;

% Pro-Saccade feature detectors.
name='ProFeature';type='exc';to='EFp';inarray=(featarray==1);retinotopic=1;
MeanInp=0.178;NoiseLevel=0;t_on=tvisual_on;t_trans_off=tvisual_off;t_off=tvisual_off;sustained_level=1;
add_input;

% Fixation feature detectors.
name='FixFeature';type='exc';to='EFf';inarray=(featarray==2);retinotopic=1;
MeanInp=0.178;NoiseLevel=0;t_on=tvisual_on;t_trans_off=tvisual_off;t_off=tvisual_off;sustained_level=1;
add_input;

% Anti-Saccade feature detectors.
name='AntiFeature';type='exc';to='EFa';inarray=(featarray==3);retinotopic=1;
MeanInp=0.178;NoiseLevel=0;t_on=tvisual_on;t_trans_off=tvisual_off;t_off=tvisual_off;sustained_level=1;
add_input;

% Space feature detectors. Used to recognize spaces at the fovea.
name='SpaceFeature';type='exc';to='EFspace';inarray=(inarray1(11)==0);retinotopic=0;
MeanInp=0.178;NoiseLevel=0;t_on=tvisual_on;t_trans_off=tvisual_off;t_off=tvisual_off;sustained_level=1;
add_input;


if select_task==7
% Anti-Saccade feature during delay
name='AntiDuringDelay';type='exc';to='EFa';inarray=inputs.external{3}.inarray;retinotopic=1;
MeanInp=0.178;NoiseLevel=0;t_on=tvisual_on+stimdelay;t_trans_off=tvisual_off+stimdelay;t_off=tvisual_off+stimdelay;sustained_level=1;
add_input;

% Visual Input during delay. This input is not really needed and could be 
% left away, if the anti-saccade input is consiedered to be not directly
% visual, or as salient as the first visual input.
name='VisualDuringDelay';type='exc';to='E4';inarray=inarray1;retinotopic=1;
MeanInp=0.056;NoiseLevel=1;t_on=tvisual_on+stimdelay;t_trans_off=tvisual_off+stimdelay;t_off=tvisual_off+stimdelay;sustained_level=0.5;
add_input;
end
