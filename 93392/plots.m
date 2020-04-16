% PLOTS.M
%
% Infer full compartmental model given only access to the voltage in the
% compartments. This code is released in conjunction with the paper 
%
%	Huys QJM, Ahrens M and Paninski L (2006): Efficient estimation of
%	detailed single-neurone models
%
% and can be downloaded from 
%
%	http://www.gatsby.ucl.ac.uk/~qhuys/code.html
%
% This scripte generates a plot of the voltage trace and the parameters
% inferred. It is called by MAIN.M. It asks whether to generate an inset to
% visualise the amount of noise in the data. If the answer is yes, you get to
% select an area in the rightmost panel to enlarge. 
%
% Copyright Quentin Huys 2006

figure(1);
clf;
subplot(2,3,1); hold on;
		ind = 1:3:nch;
		plot(atrue(ind),a(ind),'k.','markersize',10);%axis tight
		fplot('x',[min(atrue(ind)) max(atrue(ind))],'k');
		hold off
		set(gca,'fontsize',18,'xlim',[min(atrue(ind)) max(atrue(ind))],'ylim',[min(a(ind)) max(a(ind))])
		xlabel('true g_{Na}');ylabel('est g_{Na}')
		box on;
subplot(2,3,2); hold on;
		ind = 2:3:nch;
		plot(atrue(ind),a(ind),'k.','markersize',10);%axis tight
		fplot('x',[min(atrue(ind)) max(atrue(ind))],'k');
		hold off
		set(gca,'fontsize',18)
		set(gca,'fontsize',18,'xlim',[min(atrue(ind)) max(atrue(ind))],'ylim',[min(a(ind)) max(a(ind))])
		xlabel('true g_{K}');ylabel('est g_{K}')
		box on;
subplot(2,3,4); hold on;
		ind = 3:3:nch;
		plot(atrue(ind),a(ind),'k.','markersize',10);%axis tight
		fplot('x',[min(atrue(ind)) max(atrue(ind))],'k');
		hold off
		set(gca,'fontsize',18)
		set(gca,'fontsize',18,'xlim',[min(atrue(ind)) max(atrue(ind))],'ylim',[min(a(ind)) max(a(ind))])
		xlabel('true g_{L}');ylabel('est g_{L}')
		box on;
subplot(2,3,5); hold on;
		ind = nch+1:nch+nc-1;
		plot(atrue(ind),a(ind),'k.','markersize',10);%axis tight
		fplot('x',[min(atrue(ind)) max(atrue(ind))],'k');
		hold off
		set(gca,'fontsize',18)
		set(gca,'fontsize',18,'xlim',[min(atrue(ind)) max(atrue(ind))],'ylim',[min(a(ind)) max(a(ind))])
		xlabel('true g_{intercomp}');ylabel('est g_{intercomp}')
		box on;
subplot(2,3,[3 6]);	
		timevec = [1:T]*delta;
		if nc>=30; vind=1:30;
		else 	   vind = 1:nc;
		end
		plot(timevec,V(:,vind),'k');
		axis tight
		set(gca,'fontsize',18)
		xlabel('Time [ms]');ylabel('Voltage [mV]');
		box on


if str2num(input('Select region for voltage inset? [default no; yes=1]','s'))
	reg = getrect(gca);
	axes('position',[0.9050-.075    0.7    0.07    0.2]);box on
			tind = find(timevec>reg(1) & timevec<reg(1)+reg(3));
			plot(timevec(tind),V(tind,vind),'k');
			axis tight
			set(gca,'fontsize',14)
end

