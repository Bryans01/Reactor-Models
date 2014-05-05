%Inversion
clear all
close all
clc
 
%Initializing variables
Fw = []; Ff = []; Fs = []; Fg = []; Nw = []; Nf = []; Ns = []; Ng = [];
ww = []; wf = []; ws = []; wg = []; Ftot = []; Ntot = []; xw = []; xf = [];
xs = []; xg = []; Cg = []; Cf = []; Cw = []; r = []; X = []; 
%----------------------------------------------------------------------
%Start
 
%Chosen or assumed
delZ = 0.01;  %dm3 chosen
Tin = 45 + 273.15; %c 
T(1) = Tin; % 
%----------------------------------------------------------------------
Ltot = 3.75 %m
%----------------------------------------------------------------------
R = 0.5; %m radius of reactor
%----------------------------------------------------------------------
A = pi*R^2; %m2
V1 = A*delZ; %m3 is my differential volume
ee = 0.7; %void space assumed
 
densw = 1000; %kg/m3
densf = 1459;
denss = 1290;
densg = 1540;
densenz = 330;
menv = densenz*V1; %g assumed mass of enzyme
 
%Heat capacities
Cpg = 3360; %j/kg
Cps = 4290;
Cpw = 4182;
Cpf = 4120; 
 
Hr = -14930; %j/mol
%----------------------------------------------------------------------
Ftot(1) = 0.277; %mass flow rates kg/s
%----------------------------------------------------------------------
ww(1) = 0.37;
wg(1) = 0;
wf(1) = 0;
ws(1) = 0.63;
 
Fw(1)= ww(1)*Ftot(1);  %  kg/s
Ff(1)= wf(1)*Ftot(1); 
Fs(1)= ws(1)*Ftot(1);     
Fg(1)= wg(1)*Ftot(1); 
 
mmw = 18.0153*10^-3; %kg/mol
mmf = 180.156*10^-3;
mms = 342.296*10^-3;
mmg = 180.156*10^-3;
 
Nw(1) = Fw(1)/mmw; %molar flow rates mol/s
Nf(1) = Ff(1)/mmf;
Ns(1) = Fs(1)/mms;
Ng(1) = Fg(1)/mmg;
 
Qtot = (Fg(1)/densg) + (Fw(1)/densw) + (Ff(1)/densf) + (Fs(1)/denss); % m3/s
 
% j = 1 is the inlet
K = Ltot/delZ 
n = delZ/(Ltot);
Z = 1:K;
ZZ = 1:K+1;
V = Z*delZ*A;
VV = ZZ*delZ*A;
Vtot = Ltot*A;
vv = Qtot/A;
 
for j = 1:K
    Z(j) = j;
    L(j) =Z(j)*delZ;
    Tau = Qtot/Vtot;
 
    Ntot(j) = Nw(j) + Nf(j) + Ns(j) + Ng(j);
    xg(j) = Ng(j)/Ntot(j);
    xw(j) = Nw(j)/Ntot(j);
    xf(j) = Nf(j)/Ntot(j);
    xs(j) = Ns(j)/Ntot(j);
 
    Fw(j) = Nw(j)*mmw; %molar flow rates mol/s
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
    Cw(j) = Nw(j)/(Qtot); %mol/m3 
    Csu(j) = Ns(j)/(Qtot); %Cs is already taken by matlab 
 
%----------------------------------------------------------------------
    kk(j) = ((4.566*10^16)*exp(-10739/(T(j))))*10^-3*60;
    lamda(j) = (-0.003*T(j)+2.275);
    Deff(j) = ((0.403*T(j)-118.8)*10^-8)/10000;
    ro = 0.4*10^-3; %m 0.4 - 0.8 mm
    W(j) = kk(j)*(ro/Deff(j))^(1/2);
    eta(j) = (3/W(j))*(1/tanh(W(j))-1/W(j));
    
    kkk(j) = kk(j)*eta(j)*lamda(j)*(1-ee);  
%----------------------------------------------------------------------
    %Csu(j) = Csu(1)*exp(-kk(j).*Tau);
%----------------------------------------------------------------------
 
    Ns(j+1) = Ns(j)-kkk(j).*Csu(j)*V1; 
    Nw(j+1) = Nw(j)-kkk(j).*Csu(j)*V1; 
    Ng(j+1) = Ng(j)+kkk(j).*Csu(j)*V1; 
    Nf(j+1) = Nf(j)+kkk(j).*Csu(j)*V1; 
    X(j) = (Ns(1)-Ns(j+1))/Ns(1);
 
    %Nf(j+1) = Ns(j)*((Nf(1)/Ns(1))+X(j));
    %Nw(j+1) = Ns(j)*((Nw(1)/Ns(1))-X(j));
    %Ng(j+1) = Ns(j)*((Ng(1)/Ns(1))+X(j));
 
    Cpmix(j) = Cpg*xg(j)+Cpf*xf(j)+Cpw*xw(j)+Cps*xs(j);
    T(j+1) = T(j) + Hr*(Ns(j+1)-Ns(j))/(Cpmix(j)*Ntot(j));
    end
 
Fw(j+1) = Nw(j+1)*mmw; %molar flow rates kmol/hr
Ff(j+1) = Nf(j+1)*mmf;
Fs(j+1) = Ns(j+1)*mms;
Fg(j+1) = Ng(j+1)*mmg;
Ftot(j+1) = Fw(j+1) + Ff(j+1) + Fs(j+1) + Fg(j+1);
wg(j+1) = Fg(j+1)/Ftot(j+1);
ww(j+1) = Fw(j+1)/Ftot(j+1);
wf(j+1) = Ff(j+1)/Ftot(j+1);
ws(j+1) = Fs(j+1)/Ftot(j+1);
 
menvtot= menv*Ltot*R^2*pi/V1
%Conversion vs Z
subplot(2,2,1)
title('Conversion')
plot(V,X,'k')
 
%Mass flow of glucose,fructose,sucrose vs Z
subplot(2,2,2)
plot(VV,Fw,'--r',VV,Ff,'--y',VV,Fs,'--b',VV,Fg,'--c')
 
%Mol fraction of glucose,fructose,sucrose vs Z
subplot(2,2,3)
plot(VV,wg,'--r',VV,wf,'--y',VV,ws,'--b',VV,ww,'--c')
 
%Reaction rate vs Z
subplot(2,2,4)
plot(VV,T,'k')


