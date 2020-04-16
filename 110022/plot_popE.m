% plot the FOVEA, E4, E23, E5, E6 (solid line) specified in pop_active, 
% plus as comparison for a set of different positions (dotted) specified in
% pop_bg.
%
% created: Jakob Heinzle 01/07

figure;

pop_active=pp;
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
plot(E23_RATE(pop_active,t)','-','LineWidth',2);
title(['Layer 23E: Position ',(int2str(pop_active-11))]);
ylabel('frequency [Hz]');
xlim([0 length(t)]);

subplot(5,1,3);
plot(E4_RATE(pop_active,t)','-','LineWidth',2);
title(['Layer 4E: Position ',(int2str(pop_active-11))]);
ylabel('frequency [Hz]');
xlim([0 length(t)]);

subplot(5,1,4);
plot(E5R_RATE(pop_active,t)','--','LineWidth',2);
hold;
plot(E5B_RATE(pop_active,t)','-','LineWidth',2);
title(['Layer 5E: Position ',(int2str(pop_active-11))]);
ylabel('frequency [Hz]');
xlim([0 length(t)]);

subplot(5,1,5);
plot(E6A_RATE(pop_active,t)','-','LineWidth',2);
hold;
plot(E6S_RATE(pop_active,t)','--','LineWidth',2);
title(['Layer 6E: Position ',(int2str(pop_active-11))]);
ylabel('frequency [Hz]');
xlabel('time [ms]');
xlim([0 length(t)]);

