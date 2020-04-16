% plot E4, E23, E5, E6 and EREC for a single position (solid line) specified in pop_active, 
% plus as comparison for a set of different positions (dotted) specified in pop_bg.
% Plot the foveal representation of layer 2/3 and layer 6A the fixation.
%
% created: Jakob Heinzle 01/07

figure;

if ~exist('t_plot_start','var');
t_plot_start=1;
end
if ~exist('t_plot_end','var');
t_plot_end=length(E4_RATE(1,:));
end

t=t_plot_start:t_plot_end;

subplot(5,1,1);
plot(FOVEA,'-','LineWidth',2);
axis([0 t_plot_end -10.5 10]);
title(['FOVEA']);
ylabel('position');

subplot(5,1,2);
plot(E23_RATE(11,t)','-','LineWidth',2);
title(['Layer 23: Position 0']);
ylabel('frequency [Hz]');
xlim([0 length(t)]);

subplot(5,1,3);
plot(E6A_RATE(11,t)','-','LineWidth',2);
title(['Layer 6A: Position 0']);
ylabel('frequency [Hz]');
xlim([0 length(t)]);

subplot(5,1,4);
plot(IFIX_RATE,'-','LineWidth',2);
title(['Fixation Neurons']);
ylabel('frequency [Hz]');
xlim([0 length(t)]);