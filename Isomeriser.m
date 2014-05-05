%Isomerisation
clear all
close all
clc
 
%Initializing variables
Fw = []; Ff = []; Fs = []; Fg = []; Nw = []; Nf = []; Ns = []; Ng = [];
ww = []; wf = []; ws = []; wg = []; Ftot = []; Ntot = []; xw = []; xf = [];
xs = []; xg = []; Cg = []; Cf = []; r = []; 
%----------------------------------------------------------------------
%Start
 
%Chosen or assumed
delZ = 0.01;  %m3 
Tin = 50 + 273.15; %c  
T = Tin; % 
%----------------------------------------------------------------------
Ltot = 2; 
%----------------------------------------------------------------------
R = 0.125 ; %m radius of reactor
%----------------------------------------------------------------------
A = pi*R^2; %m2
V1 = A*delZ; %m3
ee = 0.4; %void space 
 
densw = 1000; %kg/m3
densf = 1459;
denss = 1290;
densg = 1540;
densenz = 330;
 
menv = densenz*V1; %g assumed mass of enzyme
 
%Heat capacities
Cpg = 336 %j/mol
Cps = 429
Cpw = 4182
Cpf = 412 
 
Hr = 3000 %j/mol
%----------------------------------------------------------------------
%Data from Eduardo Alberto
Vmf = (3.16*10^-4)/60; %mols/(s.g enz)
Vmr = (2.60*10^-4)/60; %mols/(s.g enz)
Kmf = 0.54*10^-3; %(mol/m3)
Kmr = 0.46*10^-3; %(mol/m3)
K = 1.04;
%----------------------------------------------------------------------
%j = 1 is the inlet
Ftot(1) = 13.60741 %mass flow rates kg/s
 
ww(1) = 0.476;
wg(1) = 0.388;
wf(1) = 0.049;
ws(1) = 1-ww(1)-wg(1)-wf(1);
 
Fw(1)= ww(1)*Ftot(1);  %  kg/s
Ff(1)= wf(1)*Ftot(1); 
Fs(1)= ws(1)*Ftot(1);     
Fg(1)= wg(1)*Ftot(1); 
 
mmw = 18.0153*10^3; %kg/kmol
mmf = 180.156*10^3;
mms = 342.296*10^3;
mmg = 180.156*10^3;
 
Nw(1) = Fw(1)/(mmw); %molar flow rates mol/s
Nf(1) = Ff(1)/(mmf);
Ns(1) = Fs(1)/(mms);
Ng(1) = Fg(1)/(mmg);
 
Qtot = (Fg(1)/densg) + (Fw(1)/densw) + (Ff(1)/densf) + (Fs(1)/denss);% m3/s
 
K = Ltot/delZ; 
n = delZ/(Ltot);
Z = 1:K;
ZZ = 1:K+1;
V = Z*delZ*A; % for conversion only
VV = ZZ*delZ*A; % for all else
 
for j = 1:K
    Ntot(j) = Nw(j) + Nf(j) + Ns(j) + Ng(j);
    xg(j) = Ng(j)/Ntot(j);
    xw(j) = Nw(j)/Ntot(j);
    xf(j) = Nf(j)/Ntot(j);
    xs(j) = Ns(j)/Ntot(j);
 
    Fw(j) = Nw(j)*mmw; %molar flow rates kmol/hr
    Ff(j) = Nf(j)*mmf;
    Fs(j) = Ns(j)*mms;
    Fg(j) = Ng(j)*mmg;
    Ftot(j) = Fw(j) + Ff(j) + Fs(j) + Fg(j);
    wg(j) = Fg(j)/Ftot(j);
    ww(j) = Fw(j)/Ftot(j);
    wf(j) = Ff(j)/Ftot(j);
    ws(j) = Fs(j)/Ftot(j);
 
    Cg(j) = Ng(j)/(Qtot); %mol/m3
    Cf(j) = Nf(j)/(Qtot); %mol/m3 
 
%----------------------------------------------------------------------
    r(j) = ((Kmr*Vmf*Cg(j)-Kmf*Vmr*Cf(j))/(Kmr*Kmf+Kmr*Cg(j)+Kmf*Cf(j)));
%----------------------------------------------------------------------
 
    Ng(j+1) = Ng(j)-r(j)*menv*(1-ee);
    Nf(j+1) = Nf(j)+r(j)*menv*(1-ee);
 
    X(j) = (Ng(1)-Ng(j))/Ng(1); 
 
    %Nf(j+1) = Ng(j)*((Nf(1)/Ng(1))+X(j));
 
    Ns(j+1) = Ns(j);        
    Nw(j+1) = Nw(j);   
 
    Cpmix(j) = Cpg*xg(j)+Cpf*xf(j)+Cpw*xw(j)+Cps*xs(j);
    T(j+1) = T(j) + Hr*(Ns(j+1)-Ns(j))/(Cpmix(j)*Ntot(j));
end
 
%j+1 = exit 
%To get the outlet values not just last segment values
Ntot(j+1) = Nw(j+1) + Nf(j+1) + Ns(j+1) + Ng(j+1);
 
xg(j+1) = Ng(j+1)/Ntot(j+1);
xw(j+1) = Nw(j+1)/Ntot(j+1);
xf(j+1) = Nf(j+1)/Ntot(j+1);
xs(j+1) = Ns(j+1)/Ntot(j+1);
 
Cg(j+1) = Ng(j+1)/(Qtot); %mol/m3
Cf(j+1) = Nf(j+1)/(Qtot); %mol/m3
 
Vtot = Ltot*A;
Tau= Qtot/(Vtot);
 
Fw(j+1) = Nw(j+1)*mmw; %molar flow rates kmol/hr
Ff(j+1) = Nf(j+1)*mmf;
Fs(j+1) = Ns(j+1)*mms;
Fg(j+1) = Ng(j+1)*mmg;
 
Ftot(j+1) = Fw(j+1) + Ff(j+1) + Fs(j+1) + Fg(j+1);
 
wg(j+1) = Fg(j+1)/Ftot(j+1);
ww(j+1) = Fw(j+1)/Ftot(j+1);
wf(j+1) = Ff(j+1)/Ftot(j+1);
ws(j+1) = Fs(j+1)/Ftot(j+1);
 
 
%Conversion vs Z
subplot(2,2,1)
title('Conversion')
plot(V,X,'k')
 
%Molar flow rate of glucose fructose and sucrose vs Z
subplot(2,2,2)
plot(VV,Fg,'--r',VV,Ff,'--b',VV,Fs,'--g')
 
%Mass flow rate of glucose fructose and sucrose vs Z
subplot(2,2,3)
plot(VV,wg,'--r',VV,wf,'--b',VV,ws,'--g')
 
%Mol fraction rate of glucose fructose and sucrose vs Z
subplot(2,2,4)
plot(VV,T,'k')
