function[Starts Goals Decoded Vec_l RateDiff Error] = DistanceCellModel(GC_cpm,DrawFig)
%% Distance cell / number line model of vector navigation with grid cells
%  Daniel Bush, UCL Institute of Cognitive Neuroscience
%  Reference: Using Grid Cells for Navigation (2015) Neuron (in press)
%  Contact: drdanielbush@gmail.com
%
%  Inputs:
%  GC_cpm   = Number of grid cells per unique phase, in each module
%  DrawFig  = Plot figure of errors (0 / 1)
%
%  Outputs:
%  Starts   = Random 2D start locations (m)
%  Goals    = Random 2D goal locations (m)
%  Decoded  = 2D translation vector decoded from grid cell activity (m)
%  Vec_l    = Length of decoded 2D translation vector (m)
%  RateDiff = Difference in firing rate between read-out cells on each axis
%             (Hz) - should be proportional to translation vector length
%  Error    = Error in decoded translation vector (m)

%  Provide some parameters for the simulation
iterations  = 1000;                                                         % How many iterations to run?
Range       = 500;                                                          % Range of distance cells (m)
GC_mps      = 20;                                                           % Unique grid cell phases on each axis, per module
GC_scales   = 0.25.*1.4.^(0:9);                                             % Grid cell scales (m)
GC_r        = 30;                                                           % Peak grid cell firing rate (Hz)
dt          = 0.1;                                                          % Time window of grid cell firing
Emax_k      = 0.01;                                                         % Emax WTA k parameter (see de Almeida et al., 2009)
Dist_step   = 0.04;                                                         % Displacement increment for distance cells (m)
Normalise   = 1;                                                            % Normalise distance cell firing rates? (0 / 1)
Distances   = Dist_step : Dist_step : Range;                                % Assign the locations along each axis that distance cells code for
N_distance  = length(Distances);                                            % Total number of distance cells
N_grid      = length(GC_scales)*GC_mps;                                     % Total number of grid cells scale / phase values per directional axis
clear scale

%  Generate grid cell to distance cell synaptic weights
Grid_Dist_w     = zeros(N_grid,N_distance);
for scale       = 1 : length(GC_scales)
    for offset  = 0 : GC_mps-1
        Grid_Dist_w((scale-1)*GC_mps + offset + 1, 1 : N_distance) = (cos((mod(Distances-(offset/GC_mps)*GC_scales(scale),GC_scales(scale))/GC_scales(scale))*2*pi)+1)/2;
    end
end
clear scale offset

%  Generate distance cell to read out cell synaptic weights
Dist_Out_w(1,:) = linspace(0,100,N_distance);
Dist_Out_w(2,:) = linspace(100,0,N_distance);

%  Assign some memory
Starts      = nan(iterations,2);                                            % Log start positions
Goals       = nan(iterations,2);                                            % Log goal positions
RateDiff    = nan(iterations,2);                                            % Log decoded vectors

%  For each iteration...
for i       = 1 : iterations
    
    % Update the user
    if mod(i,iterations/10)==0
        disp([int2str(i/iterations*100) '% complete...']);
        drawnow
    end
    
    % Randomly assign start and goal locations    
    Starts(i,:) = [rand*Range rand*Range];    
    Goals(i,:)  = [rand*Range rand*Range];
    
    % For each axis...
    for ax      = 1 : 2
        
        % Identify the mean firing rates of grid cells at the start and goal locations
        StartRates      = (1+cos((repmat(((mod(Starts(i,ax),GC_scales)./GC_scales)*GC_mps)',1,GC_mps) - (meshgrid(1:GC_mps,1:length(GC_scales))-1))/GC_mps*2*pi))/2*GC_r*dt;
        GoalRates       = (1+cos((repmat(((mod(Goals(i,ax), GC_scales)./GC_scales)*GC_mps)',1,GC_mps) - (meshgrid(1:GC_mps,1:length(GC_scales))-1))/GC_mps*2*pi))/2*GC_r*dt;
        
        % Convert that to Poisson firing in each of the grid cells encoding
        % each phase offset
        StartRates      = sum(poissrnd(repmat(StartRates,[1 1 GC_cpm])),3);
        GoalRates       = sum(poissrnd(repmat(GoalRates,[1 1 GC_cpm])),3);
        StartRates      = reshape(StartRates',length(GC_scales)*GC_mps,1);
        GoalRates       = reshape(GoalRates',length(GC_scales)*GC_mps,1);
    
        % Compute the firing rate of distance cells
        StartDist       = StartRates' * Grid_Dist_w;
        GoalDist        = GoalRates'  * Grid_Dist_w;
        clear StartRates GoalRates
        
        % Apply the E-max WTA algorithm
        StartDist       = StartDist .* (StartDist>((1-Emax_k).*max(StartDist)));
        GoalDist        = GoalDist  .* (GoalDist>((1-Emax_k).* max(GoalDist)));
        
        % Normalise the distance cell firing rates, if required
        if Normalise
            StartDist   = StartDist ./ sum(StartDist) * 10;
            GoalDist    = GoalDist  ./ sum(GoalDist)  * 10;
        end
        
        % Compute the firing rate of output cells
        RateDiff(i,ax)  = (sum(StartDist.*Dist_Out_w(2,:)) + sum(GoalDist.*Dist_Out_w(1,:)) - sum(StartDist.*Dist_Out_w(1,:)) - sum(GoalDist.*Dist_Out_w(2,:)))/100;
        clear StartDist GoalDist
    end                    
end

%  Compute the error
Vec_l       = sqrt(sum((Goals - Starts).^2,2));
Actual      = Goals - Starts;
b           = regress(RateDiff(:),Actual(:));
Decoded     = RateDiff./b;
Error       = sqrt(sum((Decoded - Actual).^2,2));

%  Plot the firing rate difference against true vector, if required
if DrawFig
    figure
    subplot(2,2,1)
    scatter(Actual(:),RateDiff(:),'k.')
    set(gca,'FontSize',14)
    xlabel('True Translation Vector (m)','FontSize',14)
    ylabel('Read-Out Firing Rate Difference (Hz)','FontSize',14)    
    axis square
    
    subplot(2,2,2)
    temp = histc(Error,linspace(0,ceil(max(Error)*100)/100,100)) ./ iterations;    
    bar(linspace(0,ceil(max(Error)*100),100),temp,'FaceColor','k','EdgeColor','k')
    set(gca,'FontSize',14)
    xlabel('Error in Decoded Translation Vector (cm)','FontSize',14)
    ylabel('Relative Frequency','FontSize',14)
    axis square
    
    subplot(2,2,3)
    scatter(Vec_l,Error*100,'k.')
    set(gca,'FontSize',14)
    xlabel('Decoded Translation Vector Length (m)','FontSize',14)
    ylabel('Error in Decoded Translation Vector (cm)','FontSize',14)
    hold on
    b2  = regress(Error*100,[Vec_l ones(size(Vec_l,1),1)]);
    plot(linspace(0,max(Vec_l),10),b2(2) + b2(1).*linspace(0,max(Vec_l),10),'r','LineWidth',3)
    hold off
    axis square
    clear b2
end
clear i ax Dist_Out_w Distances Grid_Dist_w iterations b N_grid N_distance Dist_step