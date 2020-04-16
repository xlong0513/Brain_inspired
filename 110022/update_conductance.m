% update_conductance updates the conductances in the network according
% to (Salinas, NeuralComputation,15,1439ff). 
%
% created: Jakob Heinzle 01/07

for k=1:cons.ncons
   if (sum(pops.population{cons.connection{k}.n_from}.spikes)~=0)
                  cons.connection{k}.G = calcg2(cons.connection{k}.G,cons.connection{k}.matrix,pops.population{cons.connection{k}.n_from}.spikes);
   end
   cons.connection{k}.G = cons.connection{k}.G*cons.connection{k}.tstep;
end