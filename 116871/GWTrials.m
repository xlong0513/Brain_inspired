% Carries out 3 trials each, with 6 different trainings, and saves
% the resulting raster plots of W1 in a series of files. Also counts
% the number of neuron firings in W1 due to C2 and C3 in the post-
% competitive epoch (from 155ms to 205ms), and saves the result in a
% file
PCE = [];
for training = 1:12
   seed1 = rand('state');
   seed2 = randn('state');
   save(['Training ',int2str(training),' rand seed'],'seed1');
   save(['Training ',int2str(training),' randn seed'],'seed2');
   GWConnect;
   for trial = 1:3
      seed1 = rand('state');
      seed2 = randn('state');
      save(['Training ',int2str(training),' Trial ',int2str(trial),...
         ' rand seed'],'seed1');
      save(['Training ',int2str(training),' Trial ',int2str(trial),...
         ' randn seed'],'seed2');
      GWRun;
      load('NetMovie','layer');
      clf;
      plot_ws_firings(layer,1,300);
      % Calculate firings due to C2 and C3 in post-competitive epoch (between
      % 150ms and 200ms)
      C2 = count_firings(layer,1,155,205,193,256);
      C3 = count_firings(layer,1,155,205,129,192);
      PCE = [PCE ; C2 C3];
      saveas(gcf,['Training ',int2str(training),' Trial ',int2str(trial)],'fig');
   end
end
save('PCE firings','PCE');