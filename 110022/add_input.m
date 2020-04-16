% add_input adds a new input to the network. 
% The following parameters need to be defined before.
% 'name': name of the input
% type: 'exc' (excitatory) or 'inh' (inhibitory)
% to: name of afferent population
% inarray: retinotopic inputarray.
% retinotopic: 1 if input is retinotopic, 0 else
% MeanInp: average input conductance
% NoiseLevel: ranges from 0 ot 1 and defines how much noise is injected.
% t_on: onset time of input
% t_trans_off: offset of phasic transient
% t_off: offset of input
% sustained_level: relative strength of sustained input to phasic input.
%
% created: Jakob Heinzle 01/07

n=inputs.ninputs+1;

n_to=0;
for k=1:pops.npops; 
   name_pop=pops.population{k}.name;
   if length(name_pop)==length(to)*strfind(name_pop,to)
   n_to=k;
   end
end
if ~n_to
   error('Invalid target population.');
end
nfac_tmp=pops.population{n_to}.nretpos;
if length(inarray)~=nfac_tmp
   error('inarray must have as many retinotopic positions as the target population.');
end

% general information about the external.
inputs.external{n}.name=name;
inputs.external{n}.retinotopic=retinotopic;
inputs.external{n}.type=type;
inputs.external{n}.to=to;
inputs.external{n}.inarray=inarray;
inputs.external{n}.MeanInp=MeanInp;
inputs.external{n}.NoiseLevel=NoiseLevel;
inputs.external{n}.t_on=t_on;
inputs.external{n}.t_off=t_off;
inputs.external{n}.t_trans_off=t_trans_off;
inputs.external{n}.sustained_level=sustained_level;

% auxiliary variables.
nperpop=pops.population{n_to}.poolsize;
InpH=zeros(nperpop*nfac_tmp,1);
for k=1:nfac_tmp
   InpH((k-1)*nperpop+1:k*nperpop)=inarray(k);
end
inputs.external{n}.ExtInp=MeanInp*InpH;
inputs.external{n}.NoiseExtInp=NoiseLevel*sqrt(gmaxE_ext*InpH/2);

pops.population{n_to}.input_external=[pops.population{n_to}.input_external, n];

inputs.ninputs=n;


 
 
 
