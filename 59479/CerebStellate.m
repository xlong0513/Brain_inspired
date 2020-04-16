
%%Model of Stellate cell, published: "A-Type and T-Type Currents Interact to Produce a Novel Spike Latency-Voltage Relationship in Cerebellar Stellate Cells.
% J Neurosci. 2005 Nov 23;25(47):10863-10873." , author Fernando R. Fernandez


clear;
warning off all;
C=1.5; %orig 1.5
Cb=2.5; %orig 0.025
s = 1;
dt=0.005;
t=0:dt:240;
I = zeros(s,length(t));
k=1:1:length(t);
kk=1:1:length(t);
kap=.25;
g=0;


gKmax=7;%7
gNamax =30;%orig 20
gNapmax =.0;%.25
gKv3max= 0; %orig 0.0
gleak=0.1; %orig .15
gA=16;%16
gC=0.55;%0.55

 td=0.01;


init =6;%at .2523
step =0.000;


Ek = -90;
ENa = 45;
Eleak = -70;
ECa=22;
Vc = zeros(1,length(t));
Vc(:,:) = -77;

Vcb = zeros(1,length(t));
Vcb(:,:) = -70;


hcinf=ones(s,length(t));

hc = zeros(s,length(t));

 
Ainf = ones(s,length(t));

AI =zeros(s,length(t));


h=1; 

n=0;
mc=0;
A=0;
m=0;

 I = zeros(s,length(t));

dVdt=zeros(s,length(t));
Vc = zeros(s,length(t));

Vc(:,:) = -80;
%     AI=1;
h=1; 

n=0;
mc=0;
A=0;
%----------------------------------------------------------------SOMA

for j = 1:1:s;
    count=1;



k=1:1:length(t);
count = 1;

step=0.1;
init=-.3;
k=1:1:length(t);

   I(j,k) = 1; 

for i = 1:1:length(t)-1;


         m = ((1/(1+exp(-(Vc(j,i)+35)/4)))).^(1/1);

         hinf = (1/(1+exp(-(Vc(j,i)+35)/-4)));
         tauh = ((2*232/pi)*(28/(4*(Vc(j,i)+74)^2 + 28^2)))-.15;
     
         h = (hinf-((hinf-h).*(exp(-(dt)./(tauh)))));
           
           A = ((1/(1+exp(-(Vc(j,i)+27)/8.8)))).^(1/1);

           
%      
%            
            Ainf(j,i+1) = (1/(1+exp(-(Vc(j,i)+68)/-6.6)));
           AI(j,i+1) = (Ainf(j,i)-((Ainf(j,i)-AI(j,i)).*(exp(-(dt)./(15)))));%15
           
            mc= ((1/(1+exp(-(Vc(j,i)+60)/3.0)))).^(1/1);

             hcinf(j,i+1) = (1/(1+exp(-(Vc(j,i)+78.0)/-3.75)));
%           
            hc(j,i+1) = (hcinf(j,i)-((hcinf(j,i)-hc(j,i)).*(exp(-(dt)./(15)))));%15
       
% 
        ninf = (1/(1+exp(-(Vc(j,i)+35)/4))).^(1);

       n = (ninf-((ninf-n).*(exp(-(dt)./( .50)))));
            


       
   VA   = ((I(j,i) - gNamax*m^1*h*(Vc(j,i) - ENa) - gC*mc^1.*hc(j,i)*(Vc(j,i) - ECa)- gA*A^1.*AI(j,i)*(Vc(j,i) - Ek)  - gKmax*n.^1*(Vc(j,i) - Ek) - gleak*(Vc(j,i) - Eleak))/C)*dt;
   VB    = ((I(j,i) - gNamax*m^1*h*((Vc(j,i)+ VA/2) - ENa) - gC*mc^1.*hc(j,i)*((Vc(j,i)+ VA/2) - ECa)- gA*A^1.*AI(j,i)*((Vc(j,i)+ VA/2) - Ek) - gKmax*n.^1*((Vc(j,i) + VA/2) - Ek) - gleak*((Vc(j,i)+VA/2) - Eleak))/C)*dt;
   VC    = ((I(j,i) - gNamax*m^1*h*((Vc(j,i)+ VB/2) - ENa) - gC*mc^1.*hc(j,i)*((Vc(j,i)+ VB/2) - ECa)- gA*A^1.*AI(j,i)*((Vc(j,i)+ VB/2) - Ek) - gKmax*n.^1*((Vc(j,i) + VB/2) - Ek) - gleak*((Vc(j,i)+VB/2) - Eleak))/C)*dt;
   VD    = ((I(j,i) - gNamax*m^1*h*((Vc(j,i)+ VC) - ENa) - gC*mc^1.*hc(j,i)*((Vc(j,i)+ VC) - ECa)- gA*A^1.*AI(j,i)*((Vc(j,i)+ VC) - Ek)- gKmax*n.^1*((Vc(j,i) + VC) - Ek) - gleak*((Vc(j,i)+VC) - Eleak))/C)*dt;
   Vc(j,i+1) = Vc(j,i) + (VA + 2*VB + 2*VC + VD)./6; 
        


 
        end;

end;



figure(1);
plot(t,Vc,'b')
hold on;


