clear all, 
close all, 
clc

Ts = 1;
z =tf('z',Ts);

%% Modelo nominal
d = 1;
Gnz = [0.7/(z-1), 0.1/(z-0.95)*z^-2;
       0.3/(z-0.9)*z^-2, 0.4/(z-0.95)];
   
Lnz = [z^-3, 0 ;
       0 , z^-2];

P = zpk(Lnz*Gnz)


%% calculo do filtro para a saída 1
z1 = 0.85; % polo do filtro

%%% sistema de equações para se obter o filtro
As1 = [0.95^2, 0.95, 1
       -2*z1, -1-z1, -2;
       1 1 1];
As1 = [0.95^2, 0.95, 1
       2, 1, 0;
       1 1 1];
bs1 = [(0.95)^3*(0.95-z1)^2;
       3*(1-z1)^2+2*(1-z1);
       (1-z1)^2];
       
X = As1\bs1
a1  = X(1)
b1 = X(2)
c1 = X(3)


Fe1 = zpk((a1*z^2 + b1*z+c1)/(z-z1)^2) % filtro da saída 1

zpk(1-z^-3*Fe1) %%% testando as condições

%% calculo do filtro para a saída 2
z2 = 0.85; % polo do filtro 

%%% sistema de equações para se obter o filtro
As2 = [0.95 1;
       1 1];
bs2 = [(0.95-z2)*0.95^2;
       (1-z2)];
       
X = As2\bs2
a2 = X(1)
b2 = X(2)
Fe2 = zpk((a2*z + b2)/(z-z2)) % filtro da saída 2

zpk(1-z^-2*Fe2) %%% testando as condições

%%% calculo das FTs S
S11 = minreal(Gnz(1,1)*(1-z^-3*Fe1))
S12 = minreal(Gnz(1,2)*(1-z^-3*Fe1))
S21 = minreal(Gnz(2,1)*(1-z^-2*Fe2))
S22 = minreal(Gnz(2,2)*(1-z^-2*Fe2))

