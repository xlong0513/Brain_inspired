% add_connection adds a new connection to the FEF network 
% The following parameters need to be defined before running this script:
% type: either 'exc' (for excitatory synapses) or 'inh' (for inhibitory
% synapses)
% from: name of the efferent population of neurons.
% to: name of the afferent population of neurons.
% weight: average synaptic weight of the connection (for full connectivity)
% weights are scaled, if a connection is sparsened.
% sparseness: percentage of active connections (between 0 and 1).  
% tau: time constant of the synapse
% small_matrix: matrix defining the population connectivity, size is number
% of afferent populations times number of efferent populations 
%
% created: Jakob Heinzle 01/07

n=cons.ncons+1;

% look for efferent population
n_from=0;
for k=1:pops.npops; 
   name_pop=pops.population{k}.name;
   if length(name_pop)==length(from)*strfind(name_pop,from)
   n_from=k;
   end
end
if ~n_from
error('Invalid from population.');
end

% look for afferent population
n_to=0;
for k=1:pops.npops; 
   name_pop=pops.population{k}.name;
   if length(name_pop)==length(to)*strfind(name_pop,to)
   n_to=k;
   end
end
if ~n_to
error('Invalid to population.');
end

% link connctions to populations.
if type=='exc'
pops.population{n_to}.input_exc=[pops.population{n_to}.input_exc, n];
elseif type=='inh'
 pops.population{n_to}.input_inh=[pops.population{n_to}.input_inh, n];
else
   error('Invalid synapse type');
end

size_to=pops.population{n_to}.nretpos;
poolsize_to=pops.population{n_to}.poolsize;
size_from=pops.population{n_from}.nretpos;
poolsize_from=pops.population{n_from}.poolsize;

if sum(size(small_matrix)==[size_to size_from])<1.6
    error('Matrix size for connection is incompatible with populations.');
end

% general parameters the connection.
cons.connection{n}.type=type;
cons.connection{n}.from=from;
cons.connection{n}.n_from=n_from;
cons.connection{n}.to=to;
cons.connection{n}.n_to=n_to;
cons.connection{n}.weight=weight;
cons.connection{n}.tau=tau;
cons.connection{n}.sparseness=sparseness;
cons.connection{n}.small_matrix=small_matrix;

%create the full weight matrix
M=zeros(poolsize_to*size_to,poolsize_from*size_from);
for toto=1:size_to;
   for fromfrom=1:size_from;
      M((toto-1)*poolsize_to+1:toto*poolsize_to,(fromfrom-1)*poolsize_from+1:fromfrom*poolsize_from)=weight*small_matrix(toto,fromfrom);
   end
end
M=sparse(M);
M=select_sparse_connect(M,sparseness);
M=distribute_random(M,w_dist,sigma_dist);

if n_to==n_from % take out autapses.
M = M-diag(diag(M));
end
cons.connection{n}.matrix=M;

% auxiliary variables of the connection
cons.connection{n}.G=zeros(pops.population{n_to}.n_neurons,1);
cons.connection{n}.tstep=exp(-dt/cons.connection{n}.tau);

cons.ncons=n;