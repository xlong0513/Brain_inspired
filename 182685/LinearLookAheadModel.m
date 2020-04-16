function[Starts Goals Decoded Vec_l Error LLA_Time] = LinearLookAheadModel(GC_cpm,DrawFig)
%% Linear look ahead model of vector navigation with grid cells
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
%  Error    = Error in decoded translation vector (m)
%  LLA_Time = Time taken to complete linear look-ahead (s)

%  Provide some parameters for the simulation
iterations  = 10;                                                         % How many iterations to run?
Range       = 500;                                                          % Size of environment (m)
GC_mps      = 20;                                                           % Unique grid cell phases on each axis, per module
GC_scales   = 0.25.*1.4.^(0:9);                                             % Grid cell scales (m)
GC_r        = 30;                                                           % Peak grid cell firing rate (Hz)
Emax_k      = 0.01;                                                         % Place cell WTA parameter
dt          = 0.005;                                                        % Place cell synaptic integration timestep (s)
SWR_speed   = 8;                                                            % Speed of linear look ahead sweep (m/s)

%  Compute the linear look ahead spatial resolution and grid module scales
Dist_step   = SWR_speed * dt;                                               % Displacement increment for place cells (m)
Distances   = 0 : Dist_step : Range;                                        % Assign place cell distance coding (m)
N_place     = length(Distances);                                            % Total number of place cells
N_grid      = length(GC_scales)*GC_mps;                                     % Total number of grid cell phase offsets
clear Dist_step

%  Generate synaptic weight matrices
Grid_Place_w        = zeros(N_grid,N_place);
for scale           = 1 : length(GC_scales)
    for offset      = 1 : GC_mps
        Grid_Place_w((scale-1)*GC_mps + offset, 1 : N_place) = ((cos((mod(Distances-((offset-1)/GC_mps)*GC_scales(scale),GC_scales(scale))/GC_scales(scale))*2*pi)+1)/2);
    end
end
clear scale offset

%  Assign some memory
Starts      = nan(iterations,2);                                            % Log of start positions
Goals       = nan(iterations,2);                                            % Log of goal positions
Decoded     = nan(iterations,2);                                            % Log of active vector cells on each axis
Error       = nan(iterations,1);                                            % Log of distance error for each computed vector
Vec_l       = nan(iterations,1);                                            % Log of true vector lengths
LLA_Time    = nan(iterations,2);                                            % Log of time taken for each linear look ahead event

%  Then, for each iteration...
for i       = 1 : iterations
    
    % Update the user
    if mod(i,iterations/10)==0
        disp([int2str(i/iterations*100) '% complete...']);
        drawnow
    end
    
    % Randomly assign start and goal locations and identify the place cells
    % encoding the goal location on each axis
    Starts(i,:) = [Range*rand Range*rand];
    Goals(i,:)  = [Range*rand Range*rand];
    Goal_ind    = [find(abs(Distances-Goals(i,1)) == min(abs(Distances-Goals(i,1)))) ...
                   find(abs(Distances-Goals(i,2)) == min(abs(Distances-Goals(i,2))))];
    
    % For each axis...
    for ax      = 1 : 2
        
        % Compute the phase of each grid cell at the starting location
        Phase   = (repmat(((mod(Starts(i,ax),GC_scales)./GC_scales)*GC_mps)',1,GC_mps) - (meshgrid(1:GC_mps,1:length(GC_scales))-1))/GC_mps*2*pi;
        Phase   = repmat(Phase,[1 1 2]);
        
        % Then run the dynamics
        finished    = [0 0];
        found       = 0;
        dirs        = [-1 1];
        t           = 1;
        while sum(finished)<2 && found==0
            
            % For each direction along the axis...
            for dir = 1 : length(dirs)
                
                % If the linear look ahead activity has not been terminated
                if finished(dir)==0
                    
                    % Compute the firing rate of place cells in each direction
                    Rates   = reshape((1+cos(Phase(:,:,dir)'))/2,length(GC_scales)*GC_mps,1);   % Compute the grid cell firing rate function
                    Rates   = sum(poissrnd(repmat(Rates*GC_r*dt,[1 1 GC_cpm])),3);              % Convert to Poisson spikes
                    Rates   = Rates' * Grid_Place_w;                                            % Convert to place cell firing rates
                    Rates   = Rates .* (Rates >= (1-Emax_k)*max(Rates));                        % Implement the WTA algorithm
                    
                    % Check the firing rate in the goal place cells and
                    % those at the end of the place cell output axis
                    if Rates(Goal_ind(ax))>0
                        Decoded(i,ax)   = dirs(dir)*(t-1)*dt*SWR_speed;                         % Record the length of the linear look ahead thus far
                        LLA_Time(i,ax)  = t*dt;
                        found           = 1;
                    elseif Rates(1) > 0 || Rates(end) > 0                                       % If linear look ahead has reached the range of the place cells
                        finished(dir)   = 1;                                                    % Terminate the linear look ahead in that direction on that axis
                    end
                    clear Rates
                    
                    if t    == length(Distances)
                        finished        = [1 1];
                    end
                    
                    % Increment the phase of grid cell firing
                    Phase(:,:,dir)      = Phase(:,:,dir) + repmat((dirs(dir)*(SWR_speed*dt)./GC_scales*2*pi)',1,GC_mps);                    
                end
            end
            t       = t + 1;
        end
        clear dirs dir t finished found Phase
    end
    clear ax Goal_ind
    
    % Compute the true vector length and error
    Error(i,1)  = sqrt(sum(((Goals(i,:) - Starts(i,:)) - Decoded(i,:)).^2,2));
    Vec_l(i,1)  = sqrt(sum((Goals(i,:) - Starts(i,:)).^2,2));
    
end

%  Plot vector length v error data, if required
if DrawFig
    figure
    subplot(2,2,1)
    temp = histc(Error,linspace(0,ceil(max(Error)*100)/100,100)) ./ iterations;    
    bar(linspace(0,ceil(max(Error)*100),100),temp,'FaceColor','k','EdgeColor','k')
    set(gca,'FontSize',14)
    xlabel('Error in Decoded Translation Vector (cm)','FontSize',14)
    ylabel('Relative Frequency','FontSize',14)
    axis square
    
    subplot(2,2,2)
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
    
    subplot(2,2,3)
    scatter([abs(Goals(:,1)-Starts(:,1)) ; abs(Goals(:,2)-Starts(:,2))],[LLA_Time(:,1) ; LLA_Time(:,2)],'k.')
    set(gca,'FontSize',14)
    xlabel('Decoded Translation Vector Length (m)','FontSize',14)
    ylabel('Time Taken to Decode Vector (s)','FontSize',14)
    axis square    
end
clear i Distances Grid_Place_w N_grid N_place dt iterations DrawFig