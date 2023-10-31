% clear all, 
% close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
Ts = 1; % perído de amostragem em segundos

%% carrega ponto de operação
pontoOperacao
hmax = 5; % nível maximo para normalização

%% carrega modelo
load('modeloDegrau')

%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 30; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 10; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

Gcoef = Gstep; % dados do modelo identificado
Nss=size(Gcoef,1); % horizonte de modelo

Gqcoef = [0;Gstepq];

umax = 0.54; % valor de controle máximo
umin = 0; % valor de controle mínimo
dumax = 0.05; % incremento de controle máximo
dumin = -dumax; % incremento de controle mínimo


%% montando as matrizes do DMC

G = zeros(N2,Nu);
G(:,1) = Gcoef(1:N2,1);


for i=2:Nu
    G(i:end,i) = G(1:end-(i-1),1);    
end

G = G(N1:end,:);

Qy = delta*eye(N2-N1+1);
Qu = lambda*eye(Nu);

Kdmc = inv(G'*Qy*G+Qu)*G'*Qy

Kdmc1 = Kdmc(1,:);

%% matrizes para o caso com restrições
Hqp = 2*(G'*Qy*G+Qu);
fqp1 = -2*G'*Qy; 

LB = repelem(dumin,Nu')';
UB = repelem(dumax,Nu')';
Rbar = tril(ones(Nu));
Rbar = [Rbar;-Rbar];

%% inicialização vetores
nin = Nss+1;
nit = 150 + nin; % número de iterações da simulação

entradas = ae0+zeros(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = h0+zeros(nit,1); % vetor com as saídas do sistema

perts = as0+zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+75:end) = as0*1.5;

refs = h0+zeros(nit,1); % vetor de referências
refs(nin+10:end) = h0*1.1;

Ylivre = ones(Nss,1)*(h0/hmax);

%% simulação sem filtro de referência
for k = nin:nit
    %% modelo processo, não mexer
    saidas(k) = simModelo(entradas(k-1),perts(k-1),saidas(k-1));
    
    %% -- Controlador DMC 
    dq = perts(k)-perts(k-1);
    Ylivre = Ylivre + Gcoef(1:Nss)*du(k-1) + Gqcoef(1:Nss)*dq;
    
    eta = (saidas(k)/hmax)-Ylivre(1); %%% erro de predição
    
    Ylivre(1:end-1) = Ylivre(2:end); %%% deslocamento do vetor
    
    %%% resposta livre
    f = Ylivre(N1:N2,1)+eta;
    
    %%% referências
    R = ones(N2-N1+1,1)*(refs(k)/hmax);
    
    
    %% Resolve o problema de otimização
    fqp = fqp1*(R-f);
    
    rbar = [repelem(umax-entradas(k-1),Nu)';
             repelem(entradas(k-1)-umin,Nu)'];
    X = quadprog(Hqp,fqp,Rbar,rbar,[],[],LB,UB);
    du(k) = X(1);
    entradas(k) = entradas(k-1)+du(k);
    
end

%% plots
t = ((nin:nit)-nin)*Ts;
vx = nin:nit;

cores = gray(3);
cores = cores(1:end-1,:);


hf = figure
h=subplot(3,1,1)
plot(t,saidas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx),'--','LineWidth',tamlinha,'Color',cores(2,:))
% ylim([0 1.5])
% h.YTick = [0 0.5 1 1.5];
hl = legend('DMC','Referência','Location','NorthEast')
ylabel(['Controlada',newline,'(m)'],'FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(3,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
% h.YTick = [-2 -1 0 1 2]
% ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on

set(h, 'FontSize', tamletra);

h = subplot(3,1,3)
plot(t,du(vx),'LineWidth',tamlinha,'Color',cores(1,:))
% h.YTick = [-2 -1 0 1 2]
% ylim([-2.5 2])
grid on
ylabel('\Delta u','FontSize', tamletra)

xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);

hf.Position = tamfigura;
hl.Position = [0.6952 0.6683 0.2054 0.1242];



