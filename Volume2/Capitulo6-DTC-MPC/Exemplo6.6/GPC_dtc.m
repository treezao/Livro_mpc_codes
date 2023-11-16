clear all, 
close all, 
clc

%% obtenção do modelo

Ts = 0.1; % período de amostragem
z =tf('z',Ts);
s = tf('s');

Ke = 2;
tau = 1;
G = Ke/(tau*s+1)^2
L1 = 5*Ts;

P = G*exp(-L1*s);

Gnz = c2d(G,Ts,'zoh')
d = L1/Ts;
Lnz = z^-d;


%% ajuste GPC
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

m = 1; % número de entradas da planta - manipuladas
n = 1; % número de saidas da planta - controladas

Nu = 2; %horizontes de controle - dimensão 1 x m
N1 = d+1; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = d+3; % fim dos horizontes de predição - dimensão 1 x n
Ny = N2-N1+1; %horizonte de predição - não editar

delta = 1; %ponderação nos erros - dimensão 1 x n
lambda = 1; %ponderação nas ações de controle - dimensão 1 x m       

%%% montagem das matrizes de predição do DTC-GPC

nH = zeros(1,m); % qt de dados passados necessários das entradas
nF = zeros(1,n); % qt de dados passados necessários das saídas.

[Ei,Fi] = diofantina(conv(Gnz.den{1},[1 -1]),N1-d,N2-d);
    
F = Fi;
nF = size(Fi,2);
    
% incorporação do atraso no numerador B
Btil = Gnz.num{1}(2:end); 

Gtemp = [];
Htemp = [];
for k=1:N2-d
    EjB = conv(Ei(k,1:k),Btil);

    Gtemp(k,1:k) = EjB(k:-1:1);    
    Htemp(k,:) = EjB(k+1:end);
end
G = Gtemp(:,1:Nu);
H = Htemp;

nH = size(H,2);
        
%%% matrizes de ponderação
Qy = diag(repelem(delta,1,Ny)); %montagem da matriz de ponderação do erro
Qu = diag(repelem(lambda,1,Nu)); %montagem da matriz de ponderação da ação de controle


Kgpc = (G'*Qy*G+Qu)\G'*Qy;
Kgpc1 = Kgpc(1,:);

%%

f0 = sum(Kgpc1)
beq = Kgpc1*H
aeq = Kgpc1*F

Dc = conv([1,-1],[1 beq]);
Nc = aeq;


[Ei,Fi] = diofantina(conv(Gnz.den{1},[1 -1]),1,d);

Nfe = 0;
for i=1:size(Fi,2)
    
    if(d-i+1>0)
        Nfe = Nfe + aeq(i)*Fi(d-i+1,:);
    else
        temp = zeros(1,size(Fi,2));
        temp(-(d-i+1)+1)=1;
        Nfe = Nfe + temp;
    end
end

Fref = tf(f0,Nfe,Ts,'Variable','z^-1')
zpk(Fref)
Fe = tf(Nfe,Nc,Ts,'Variable','z^-1')
zpk(Fe)

C = tf(Nc,Dc,Ts,'Variable','z^-1')
zpk(C)

