% Sollini et al. ON-OFF receptive fields in auditory cortex diverge during
% development and contribute to directional sweep selectivity. Nat. Com
% 2018
clear all
% Parameters
T = 100000;
Ninp = 20;
alphae = 0.0001;
alphai = 0.00001;
rho = 0.01;
wmax = 1;
win = 1;
Sig = 1.5;
ythres = 2.5;
L1 = 1/Ninp;
L1i = 1/Ninp/2;
wnoise = 0.005;

% Create currents
G = zeros(Ninp/2, Ninp/2);
G(5,:) = exp(-((Ninp/4)-(1:Ninp/2)).^2/(2*Sig^2));
for i = 2:Ninp/2
    G(mod(i+3,Ninp/2)+1,:) = circshift(G(mod(i+2,Ninp/2)+1,:)',1)';
end
G = [G,G*0;G*0,G];


% Initialisation of weights
we = win*rand(Ninp,1)/Ninp;
wi = 0*win*rand(Ninp,1)/Ninp;

% Create inputs before hearing onset
tau_OU = 5;
xe = zeros(Ninp/2,T);
for t = 2:T
    xe(:,t) = xe(:,t-1)*exp(-1/tau_OU)+(0.5-rand(Ninp/2,1))*(1-exp(-1/tau_OU));
end
xe = (xe-0.1).*(xe>0.1);
xe = xe/mean(mean(xe));
xbar  = [xe ;xe];

% Simulations before hearing onset
for t = 2:T
    y = we'* G*xbar(:,t)- wi'* G*xbar(:,t);
    y = ((y-ythres)>0)*(y-ythres);
    we = we + alphae*G*xbar(:,t)*y+wnoise*we.*(rand(Ninp,1)-0.5);
    we(1:Ninp/2) = we(1:Ninp/2) - sum(we(1:Ninp/2))/Ninp/2+L1;
    we(1+Ninp/2:end) = we(1+Ninp/2:end) - sum(we(1+Ninp/2:end))/Ninp/2+L1;
    we = we.*(we>0);
    we = (we-wmax).*(we<wmax)+wmax;
    wi = wi + alphai*G*xbar(:,t)*(y-rho);
    wi(1:Ninp/2) = wi(1:Ninp/2) - sum(wi(1:Ninp/2))/Ninp/2+L1i;
    wi(1+Ninp/2:end) = wi(1+Ninp/2:end) - sum(wi(1+Ninp/2:end))/Ninp/2+L1i;
    wi = wi.*(wi>0);
end


% Create Inputs after hearing onset - sound evoked
rate     = 50;
tau_onff = 10;
xe = zeros(Ninp/2,T);
xbar = zeros(Ninp,T);
freq_index = floor(Ninp/2*rand)+1;
for t = 2:T
    bswitch = rand<(1/rate);
    if bswitch ==1;if xe(freq_index,t-1)==0; freq_index = floor(Ninp/2*rand)+1; end; end
    xbar(freq_index,t)=xbar(freq_index,t-1)*exp(-1/tau_onff)+ bswitch.*(1-xe(freq_index,t-1))*(1-exp(-1/tau_onff));
    xbar(freq_index+(Ninp/2),t)= xbar(freq_index+(Ninp/2),t-1)*exp(-1/tau_onff)  + bswitch.*xe(freq_index,t-1)*(1-exp(-1/tau_onff));
    xbar([1:freq_index-1,freq_index+1:Ninp/2,[1:freq_index-1,freq_index+1:Ninp/2]+Ninp/2],t)= xbar([1:freq_index-1,freq_index+1:Ninp/2,[1:freq_index-1,freq_index+1:Ninp/2]+Ninp/2],t-1)*exp(-1/tau_onff);
    xe(freq_index,t) = xe(freq_index,t-1)+ bswitch.*(1-2*xe(freq_index,t-1));
end
xbar  = xbar/max(max(xbar))*40;

% Simulations after hearing onset - sound evoked
for t = 2:T
    y = we'* G*xbar(:,t)- wi'* G*xbar(:,t);
    y = ((y-ythres)>0)*(y-ythres);
    we = we + alphae*G*xbar(:,t)*y+wnoise*we.*(rand(Ninp,1)-0.5);
    we(1:Ninp/2) = we(1:Ninp/2) - sum(we(1:Ninp/2))/Ninp/2+L1;
    we(1+Ninp/2:end) = we(1+Ninp/2:end) - sum(we(1+Ninp/2:end))/Ninp/2+L1;
    we = we.*(we>0);
    we = (we-wmax).*(we<wmax)+wmax;
    wi = wi + alphai*G*xbar(:,t)*(y-rho);
    wi(1:Ninp/2) = wi(1:Ninp/2) - sum(wi(1:Ninp/2))/Ninp/2+L1i;
    wi(1+Ninp/2:end) = wi(1+Ninp/2:end) - sum(wi(1+Ninp/2:end))/Ninp/2+L1i;
    wi = wi.*(wi>0);
end

figure; hold on
plot(G(1:Ninp/2,1:Ninp/2)*we(1:Ninp/2), 'linewidth', 3)
plot(G(1:Ninp/2,1:Ninp/2)*we(1+Ninp/2:end), 'linewidth', 3)
xlabel('frequency index','fontsize',20)
xlim([1 10])
ylabel('e current [a.u.]','fontsize',20)
legend('On', 'Off')
set(gca,'fontsize',20);
set(gca,'fontsize',20);