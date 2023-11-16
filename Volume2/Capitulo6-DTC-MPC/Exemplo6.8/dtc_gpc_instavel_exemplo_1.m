clear all, close all, clc
%% Exemplo instável DTC-GPC
Ts = 1;
z = tf('z',Ts)

z0 = 1.1;
d = 5;
Gz = 0.2/(z-z0)
Pz = Gz*z^-d

%%% resolvendo o sistema de equações:
zf = 0.9;

Asis = [1 1;
        z0 1]
bsis = [(1-zf)*1;
        (z0-zf)*z0^d]
    
X = Asis\bsis;
af0 = X(1)
af1 = X(2);

Fe = (af0*z+af1)/(z-zf)
zpk(Fe)

%%% obtendo o polinômio S
1-z0^-d*(af0*z0+af1)/(z0-zf)


S = Gz*(1-z^-d*Fe)
Sx = minreal(S)

Sx.num{1}./Sx.num{1}(2)

%%%
