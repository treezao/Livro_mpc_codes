clear all, 
close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%

Ts = 0.15;
z =tf('z',Ts);
s = tf('s');

%% caso 1
Ke = 1;
tau = 1.5;
G = Ke/(tau*s+1)
L1 = 10*Ts;

P = G*exp(-L1*s);

%%% ajuste PS - controlador primário
Czps = 3.70*(z-0.77)/(z-1);
Fzps = 0.23/(z-0.77);

Gnz = c2d(G,Ts,'zoh')
d = L1/Ts;
Lnz = z^-d;


%% ajuste GPC
%%% PARAMETROS DE AJUSTE DO CONTROLADOR

m = 1; % número de entradas da planta - manipuladas
n = 1; % número de saidas da planta - controladas

Nu = 5; %horizontes de controle - dimensão 1 x m
N1 = d+1; % inicio dos horizontes de predição - incluir atraso - dimensão 1 x n
N2 = d+20; % fim dos horizontes de predição - dimensão 1 x n
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

%% dados dos critérios de robustez
w = logspace(-1,log10(2*pi*1/2/Ts),200);
[mag,pha] = bode([Fe,tf(1,1,Ts)],w);

mag1 = reshape(mag(1,1,:),200,1);
mag2 = reshape(mag(1,2,:),200,1);

%%
cores = gray(4);
cores = cores(1:end-1,:);

hf = figure
h= subplot(1,1,1)
semilogx(w,20*log10(mag1),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
semilogx(w,20*log10(mag2),'--','LineWidth',tamlinha,'Color',cores(2,:))
semilogx(2*pi/2/Ts*[1 1],[-500 500],'LineWidth',tamlinha,'Color',cores(1,:))
grid on

hl = legend('F_e(z) - GPC','F_e(z) - PS','Location','NorthWest')
ylabel('Magnitude (dB)','FontSize', tamletra)

% xlim([10^-2 3*10^1])
ylim([-1 6])

set(h, 'FontSize', tamletra);


xlabel('Frequência (rad/min)','FontSize', tamletra)


set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
% hl.Position = [0.6845 0.5863 0.2054 0.1806]; 
% print('trocador_PS_GPC_dtc','-depsc')
