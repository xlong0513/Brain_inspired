function[Starts Goals Decoded Vec_l FirstError LastError Steps] = VectorCellModel(GC_cpm,DrawFig)
%% Rate-coded vector cell model of vector navigation with grid cells
%  Daniel Bush, UCL Institute of Cognitive Neuroscience
%  Reference: Using Grid Cells for Navigation (2015) Neuron (in press)
%  Contact: drdanielbush@gmail.com
%
%  Inputs:
%  GC_cpm       = Number of grid cells per unique phase, in each module
%  DrawFig      = Plot figure of errors (0 / 1)
%
%  Outputs:
%  Starts       = Random 2D start locations (m)
%  Goals        = Random 2D goal locations (m)
%  Decoded      = 2D translation vector decoded from grid cell activity (m)
%  Vec_l        = Length of decoded 2D translation vector (m)
%  FirstError   = Error in first decoded translation vector (m)
%  LastError    = Error in final decoded translation vector (m)
%  Steps        = No. of iterative steps used to compute translation vector

%  Provide some parameters for the simulation
iterations  = 1000;                                                         % How many iterations to run?
Range       = 500;                                                          % Range (m)
GC_mps      = 20;                                                           % Unique grid cell phases on each axis, per module
GC_scales   = 0.25.*1.4.^(0:9);                                             % Grid cell scales (m)
GC_r        = 30;                                                           % Peak grid cell firing rate (Hz)
dt          = 0.1;                                                          % Time window of grid cell firing (s)
Emax_k      = 0.01;                                                         % E-max winner-take-all parameter (see de Almeida et al., 2009)
N_vec_f     = 12500;                                                        % Number of fine-grained vector cell dendrites
N_vec_c     = 1250;                                                         % Number of coarse-grained vector cells

%  Assign the vectors encoded by the vector cell dendrites and cells
Vectors_f   = linspace(0,Range,N_vec_f+1);                                  % Assign the linear range of vectors encoded by the fine-grained vector cells (dendrites)
Vectors_c   = 0;                                                            % Assign the psuedo-exponential range of vectors encoded by the coarse-grained vector cells
for scale   = 1 : length(GC_scales)
    Vectors_c   = [Vectors_c Vectors_c(end)+linspace(GC_scales(scale).*(Range/(sum(GC_scales)*100)),GC_scales(scale).*Range/sum(GC_scales),round(N_vec_c/length(GC_scales)))];
end
clear scale

%  Generate synaptic weight matrices (grid cells to vector cell dendrites /
%  vector cell dendrites to cells)
Grid_Vec_w  = zeros(GC_mps,GC_mps,GC_mps);                                  % Multiplicative grid cell output synapses
for offset  = 1 : GC_mps
    for vec = 1 : GC_mps
        syn = offset + vec - 1;
        syn(syn>GC_mps)             = syn-GC_mps;
        Grid_Vec_w(offset,syn,vec)  = 5e-4;
        clear syn
    end
end
clear offset vec
Vec_Vec_w       = zeros(GC_mps,N_vec_f+1,length(GC_scales));                % Grid cell to fine-grained vector cell dendrite synapses
for scale       = 1 : length(GC_scales)
    for offset  = 1 : GC_mps
        Vec_Vec_w(offset,:,scale)   = (cos((mod(Vectors_f-((offset-1)/GC_mps)*GC_scales(scale),GC_scales(scale))/GC_scales(scale))*2*pi)+1)/2*5e-4;
    end
end
clear scale offset
Vec_fc_w        = zeros(N_vec_f,N_vec_c);                                   % Fine-grained vector cell dendrite to coarse-grained vector cell synapses
for c           = 1 : length(Vectors_f)
    [i ind]     = min(abs(Vectors_f(c)-Vectors_c));
    Vec_fc_w(c,ind) = 1;
    clear i ind
end
clear c

%  Assign some memory for the output
Starts      = nan(iterations,2);                                            % Log of start positions
Goals       = nan(iterations,2);                                            % Log of goal positions
Decoded     = nan(iterations,2,5);                                          % Log of active vector cells on each axis
Vec_l       = nan(iterations,1);                                            % Log of true vector lengths
FirstError  = nan(iterations,1);                                            % Log of first distance error for each computed vector
LastError   = nan(iterations,1);                                            % Log of last distance error for each computed vector
Steps       = nan(iterations,1);                                            % Log of steps taken

%  For each iteration...
for i       = 1 : iterations
    
    % Update the user
    if mod(i,iterations/10)==0
        disp([int2str(i/iterations*100) '% complete...']);
        drawnow
    end
    
    % Randomly assign start and goal locations
    Starts(i,:)     = [Range*rand Range*rand];
    Goals(i,:)      = [Range*rand Range*rand];
    
    % Start the iterative vector navigation process
    step            = 1;
    stop            = false;
    while ~stop
        
        % For each axis...
        for ax          = 1 : 2
            
            % ...assign memory for the output
            vector_out  = zeros(N_vec_f+1,2);
            if step     == 1
                start   = Starts(i,:);
            end
            
            % Then, for each scale...
            for scale   = 1 : length(GC_scales)
                
                % Generate the mean firing rate of all cells in that grid module on that axis
                StartRates      = (1+cos(((mod(start(1,ax),GC_scales(scale))/GC_scales(scale))*GC_mps-(1:GC_mps)')/GC_mps*2*pi))/2*GC_r*dt;
                GoalRates       = (1+cos(((mod(Goals(i,ax),GC_scales(scale))/GC_scales(scale))*GC_mps-(1:GC_mps)')/GC_mps*2*pi))/2*GC_r*dt;
                
                % Convert that to Poisson firing in each of the grid cells encoding each phase offset
                StartRates      = sum(poissrnd(repmat(StartRates,[1 GC_cpm])),2);
                GoalRates       = sum(poissrnd(repmat(GoalRates,[1 GC_cpm])),2);
                
                % Compute the output of each grid cell module through the multiplicative synapses
                multsyn_dir1    = sum(squeeze(sum(repmat(StartRates * GoalRates', [1 1 GC_mps]).*Grid_Vec_w)))';
                multsyn_dir2    = sum(squeeze(sum(repmat(GoalRates  * StartRates',[1 1 GC_mps]).*Grid_Vec_w)))';
                
                % Compute the input to each fine grained vector cell dendrite from that grid cell module
                input_dir1      = (multsyn_dir1' * Vec_Vec_w(:,:,scale))';
                input_dir2      = (multsyn_dir2' * Vec_Vec_w(:,:,scale))';
                
                % Store the overall output of vector cell dendrites
                vector_out(:,1) = vector_out(:,1)   + input_dir1;
                vector_out(:,2) = vector_out(:,2)   + input_dir2;
                clear StartRates GoalRates multsyn_dir1 multsyn_dir2 input_dir1 input_dir2
            end
            
            % Implement the WTA algorithm and decode the vector as the weighted
            % mean of all vector cells firing in each array
            vector_out          = double(vector_out>((1-Emax_k)*max(vector_out(:))));
            Decoded(i,ax,step)  = nanmean([Vectors_c((vector_out(:,1)'*Vec_fc_w)>0) -Vectors_c((vector_out(:,2)'*Vec_fc_w)>0)]);
            
        end
        clear ax scale
        
        % Then, unless the last position is within one metre of the
        % goal, move to the new start position and take another step
        if (sqrt(sum((start + 0.8 .*Decoded(i,:,step) - Goals(i,:)).^2)) > 1) && (step <= 5)
            start           = start + 0.8 .* Decoded(i,:,step);
            step            = step + 1;
        else
            stop            = true;
        end
        clear vector_out
        
    end
    
    % Compute the true vector length, first step error and total number of steps
    FirstError(i,1) = sqrt(sum(((Goals(i,:) - Starts(i,:)) - Decoded(i,:,1)).^2,2));
    LastError(i,1)  = sqrt(sum((Goals(i,:) - start - Decoded(i,:,step)).^2,2));
    Vec_l(i,1)      = sqrt(sum((Goals(i,:) - Starts(i,:)).^2,2));
    Steps(i,1)      = step;
    clear step start stop
    
end

%  Plot vector length v error data, if required
if DrawFig
    figure
    subplot(2,2,1)
    scatter(1:length(Vectors_c),Vectors_c,'k.')
    set(gca,'FontSize',14)
    xlabel('Vector Cell Index','FontSize',14)
    ylabel('Encoded Vector (m)','FontSize',14)
    axis square
    
    subplot(2,2,2)
    temp = histc(LastError,linspace(0,ceil(max(LastError)*100)/100,100)) ./ iterations;    
    bar(linspace(0,ceil(max(LastError)*100),100),temp,'FaceColor','k','EdgeColor','k')
    set(gca,'FontSize',14)
    xlabel('Error in Decoded Translation Vector (cm)','FontSize',14)
    ylabel('Relative Frequency','FontSize',14)
    axis square
    
    subplot(2,2,3)
    scatter(Vec_l,FirstError*100,'k.')
    set(gca,'FontSize',14)
    xlabel('True Translational Vector Length (m)','FontSize',14)
    ylabel('First Decoded Translational Vector Error (cm)','FontSize',14)    
    b2 = regress(FirstError*100,[Vec_l ones(size(Vec_l,1),1)]);
    hold on    
    plot(linspace(0,max(Vec_l),10),b2(2) + b2(1).*linspace(0,max(Vec_l),10),'r','LineWidth',3)
    hold off
    axis square
    clear b2
end
clear i Grid_Vec_w Vec_Vec_w Emax_k n_iters N_vec_c N_vec_f Norm Range iterations Vec_fc_w Vectors_f