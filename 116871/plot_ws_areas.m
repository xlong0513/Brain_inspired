function plot_ws_areas(layer,i,j,tmax)
% Produces side by side raster plots of workspace areas i and j
subplot(1,2,1)
plot_ws_firings(layer,i,tmax);
subplot(1,2,2)
plot_ws_firings(layer,j,tmax);
set(gcf,'Position',[50,100,1200,400]);