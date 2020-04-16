% add_population adds a new population to the FEF network
% the following paramters need to be defined earlier:
% name: name of the population
% type: 'exc' (excitatory) or 'inh' (inhibitory)
% poolsize: size of one population in the array
% nretpos: number of retinotopic positions
% bgE: mean excitatory background input.
% bgI: mean inhibitory background input.
%
% created: Jakob Heinzle 01/07

n=pops.npops+1;

% general information about the population.
pops.population{n}.name=name;
pops.population{n}.type=type;
pops.population{n}.poolsize=poolsize;
pops.population{n}.nretpos=nretpos;
pops.population{n}.bgExc=bgE;
pops.population{n}.bgInh=bgI;
pops.population{n}.NoiseExc=sqrt(gmaxE_ext*pops.population{n}.bgExc/2);
pops.population{n}.NoiseInh=sqrt(gmaxI_ext*pops.population{n}.bgInh/2);

% auxiliary variables.
n_neurons=pops.population{n}.nretpos*pops.population{n}.poolsize;
pops.population{n}.n_neurons=n_neurons;
pops.population{n}.tspk=zeros(n_neurons,1);
pops.population{n}.Vm=5*rand(n_neurons,1);
pops.population{n}.spikes=zeros(n_neurons,1);
pops.population{n}.refrac = zeros(n_neurons,1);
pops.population{n}.inEAux = ones(n_neurons,1)*pops.population{n}.bgExc;
pops.population{n}.inIAux = ones(n_neurons,1)*pops.population{n}.bgInh;

% information added later, when connections and external inputs are added.
pops.population{n}.input_external=[]; %define, which external inputs, target the population
pops.population{n}.input_exc=[]; %define, which internal conductances, target the population
pops.population{n}.input_inh=[]; %define, which internal conductances, target the population

pops.npops=n;
