function[Starts Goals Decoded Vec_l FirstError LastError Steps] = PhaseModel(PNoise,GC_cpm,DrawFig)
%% Phase coded vector cell model of vector navigation with grid cells
%  Daniel Bush, UCL Institute of Cognitive Neuroscience
%  Reference: Using Grid Cells for Navigation (2015) Neuron (in press)
%  Contact: drdanielbush@gmail.com
%
%  Inputs:
%  PNoise       = Standard deviation of phase noise (rad)
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
dt          = 0.1;                                                          % Length of a theta cycle (s)
Emax_k      = 0.01;                                                         % E-max winner-take-all parameter (see de Almeida et al., 2009)
N_vec_f     = 12500;                                                        % Number of fine-grained vector cells (dendrites?)
N_vec_c     = 1250;                                                         % Number of coarse-grained vector cells

%  Assign the vectors encoded by the vector cells
Vectors_f   = linspace(0,Range,N_vec_f+1);                                  % Assign the linear range of vectors encoded by the fine-grained vector cells (dendrites)
Vectors_c   = 0;                                                            % Assign the psuedo-exponential range of vectors encoded by the coarse-grained vector cells
for scale   = 1 : length(GC_scales)
    Vectors_c   = [Vectors_c Vectors_c(end)+linspace(GC_scales(scale).*(Range/(sum(GC_scales)*100)),GC_scales(scale).*Range/sum(GC_scales),round(N_vec_c/length(GC_scales)))];
end
clear scale

%  Compute grid cell peak firing locations and generate the delay line connectivity
Locations   = (meshgrid(1:GC_mps,1:length(GC_scales))-1) / GC_mps .* repmat(GC_scales',1,GC_mps);
DelayLines  = nan(length(GC_scales),length(Vectors_f));
for scale   = 1 : length(GC_scales)
    DelayLines(scale,:) = mod(GC_scales(scale)/2 - Vectors_f,GC_scales(scale))/GC_scales(scale)*dt;
end
clear scale

%  Generate the the fine to coarse grained vector cell synaptic connectivity
Vec_fc_w    = zeros(N_vec_f,N_vec_c);
for c       = 1 : length(Vectors_f)
    [i ind] = min(abs(Vectors_f(c)-Vectors_c));
    Vec_fc_w(c,ind)     = 1;
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
Steps       = nan(iterations,1);                                            % Number of steps taken

%  Then, for each iteration...
for i       = 1 : iterations
    
    % Update the user
    if mod(i,iterations/10)==0
        disp([int2str(i/iterations*100) '% complete...']);
        drawnow
    end
    
    % Choose a random start and goal location
    Starts(i,:) = [Range*rand Range*rand];
    Goals(i,:)  = [Range*rand Range*rand];
    
    % Start the iterative vector navigation process
    start       = Starts(i,:);
    step        = 1;
    stop        = false;
    while ~stop
        
        % Then, for each axis...
        for ax      = 1 : 2
            
            % Compute the phase of firing in each grid cell at the current
            % location, in each direction along the axis
            StartPhases1        = mod(mod((meshgrid(1:GC_mps,1:length(GC_scales))-1) / GC_mps .* repmat(GC_scales',1,GC_mps) - start(ax), repmat(GC_scales',1,GC_mps)) ./ repmat(GC_scales',1,GC_mps) * 2*pi + pi, 2*pi);
            StartPhases2        = mod(mod(start(ax) - (meshgrid(1:GC_mps,1:length(GC_scales))-1) / GC_mps .* repmat(GC_scales',1,GC_mps), repmat(GC_scales',1,GC_mps)) ./ repmat(GC_scales',1,GC_mps) * 2*pi + pi, 2*pi);
            
            % Compute the set of grid cells across modules that fire maximally at the goal
            [Dist Cells]        = min(abs(Locations - repmat(mod(Goals(i,ax),GC_scales)',1,GC_mps)),[],2);
            Inds                = sub2ind(size(Locations),(1:length(GC_scales))',Cells);
            clear Dist Cells
            
            % Compute the firing times of that set of grid cells across modules
            FiringTimes1        = (repmat(StartPhases1(Inds),GC_cpm,1) + PNoise*randn(length(GC_scales)*GC_cpm,1))/(2*pi)*dt;
            FiringTimes2        = (repmat(StartPhases2(Inds),GC_cpm,1) + PNoise*randn(length(GC_scales)*GC_cpm,1))/(2*pi)*dt;
            clear Inds StartPhases1 StartPhases2
            
            % Identify the coincidence with which these spikes arrive at each
            % fine-grained vector cell (dendrite), in each direction
            ArrivalTimes        = [abs(sum(exp(1i.*(mod(repmat(DelayLines,GC_cpm,1) + repmat(FiringTimes1,1,N_vec_f+1),dt)/dt * 2 * pi))))./length(FiringTimes1) ; ...
                                   abs(sum(exp(1i.*(mod(repmat(DelayLines,GC_cpm,1) + repmat(FiringTimes2,1,N_vec_f+1),dt)/dt * 2 * pi))))./length(FiringTimes2)];
            clear FiringTimes1 FiringTimes2
                        
            % Implement the WTA algorithm and decode the vector as the weighted
            % mean of all vector cells firing in each array
            vector_out          = double(ArrivalTimes>((1-Emax_k)*max(ArrivalTimes(:))));
            Decoded(i,ax,step)  = nanmean([Vectors_c((vector_out(1,:)*Vec_fc_w)>0) -Vectors_c((vector_out(2,:)*Vec_fc_w)>0)]);
            clear ArrivalTimes vector_out            
            
        end
        clear ax
        
        % Then, unless the last position is within one metre of the
        % goal, move to the new start position and take another step
        if (sqrt(sum((start + 0.8 .*Decoded(i,:,step) - Goals(i,:)).^2)) > 1) && (step <= 5)
            start   = start + 0.8 .* Decoded(i,:,step);
            step    = step + 1;
        else
            stop    = true;
        end
        
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
    position    = repmat(linspace(-0.5,0.5,100),100,1);
    phase       = -position*2*pi + PNoise*randn(100,100);
    scatter(position(:),phase(:),'k.')
    set(gca,'FontSize',14)
    xlabel('Distance through the grid field (au)','FontSize',14)
    ylabel('Theta firing phase (rad)','FontSize',14)
    axis square
    clear phase position
    
    subplot(2,2,2)
    scatter(mod(Vectors_f+GC_scales(1)/2,GC_scales(1))/GC_scales(1)-0.5,DelayLines(1,:)*1000,'k.','SizeData',800)    
    set(gca,'FontSize',14)    
    xlabel('Relative Position of Goal (s_i)','FontSize',14)
    ylabel('Delay Line Length (ms)','FontSize',14)    
    axis square    
    
    subplot(2,2,3)
    temp = histc(LastError,linspace(0,ceil(max(LastError)*100)/100,100)) ./ iterations;    
    bar(linspace(0,ceil(max(LastError)*100),100),temp,'FaceColor','k','EdgeColor','k')
    set(gca,'FontSize',14)
    xlabel('Error in Decoded Translation Vector (cm)','FontSize',14)
    ylabel('Relative Frequency','FontSize',14)
    axis square
    
    subplot(2,2,4)
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
clear i iterations DelayLines Locations dt error N_vec_f N_vec_c Vec_fc_w Vectors_f