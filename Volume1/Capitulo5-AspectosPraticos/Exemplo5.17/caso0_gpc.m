% clear all, 
% close all, 
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%
Ts = 1; % período de amostragem
A = conv([1 -0.9],[1 -.7]); % denominador do modelo do processo
B = 0.5; % numerador do modelo do processo
Bq = [0 B]; % numerador do modelo da perturbação do processo

d=5; % atraso em relação à manipulada
dq = 5; % atraso em relação à perturbação

nb = size(B,2)-1; % ordem do polinômio B
na = size(A,2)-1; % ordem do polinômio A
nbq = size(Bq,2)-1; % ordem do polinômio Bq



%% parâmetros de ajuste

N1 = 6; %horizonte de predição inicial
N2 = 30; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 5; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

umax = 1; % valor de controle máximo
umin = 0; % valor de controle mínimo
dumax = 0.05; % incremento de controle máximo
dumin = -dumax; % incremento de controle mínimo


%% montagem das matrizes
Btil = conv(B,[zeros(1,d) 1]); % incorporação do atraso no numerador B


G = zeros(N,N); % matriz dinâmica G
H = zeros(N,nb+d); % matriz dos incrementos passados de controle 

[E,F] = diofantina(conv(A,[1 -1]),N1,N2); % obtenção dos polinômios Ej, Fj

for i=N1:N2
    EjB = conv(E(i-N1+1,1:i),Btil);
    
    G(i-N1+1,1:i) = EjB(i:-1:1);    
    H(i-N1+1,:) = EjB(i+1:end);

    
end
G = G(:,1:Nu);


G,F,H

%% obtenção do ganho irrestrito
Qu = lambda*eye(Nu); % matriz de ponderação dos incrementos de controle
Qy = delta*eye(N);  % matriz de ponderaão dos erros futuros

Kmpc = (G'*Qy*G+Qu)\G'*Qy;

Kmpc1 = Kmpc(1,:);

%% matrizes para o caso com restrições
Hqp = 2*(G'*Qy*G+Qu);
fqp1 = -2*G'*Qy; 

LB = repelem(dumin,Nu')';
UB = repelem(dumax,Nu')';
Rbar = tril(ones(Nu));
Rbar = [Rbar;-Rbar];

%% inicialização vetores
nin = 10;
nit = 30 + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
% perts(nin+50:end) = 0.5;

refs = 0*ones(nit+N2,1); % vetor de referências
refs(nin+15:end) = 1.5;



%% simulação sem filtro de referência
for k = nin:nit
    %% modelo processo, não mexer
    saidas(k) = -A(2:end)*saidas(k-1:-1:k-na) ...
                  +B*entradas(k-d-1:-1:k-1-nb-d) ...
                  +Bq*perts(k-dq:-1:k-nbq-dq);    
    
    %% -- Controlador GPC 
    %%% referencias
    R = refs(k+N1:k+N2);
    
    %%% cálculo da resposta livre;
    f = F*saidas(k:-1:k-na);
    
    if(~isempty(H))
        f = f+ H*du(k-1:-1:k-nb-d); % parcela dos incrementos de controle
    end
    
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
hl = legend('GPC','Referência','Location','NorthEast')
ylabel('Controlada','FontSize', tamletra)
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

% hf.Position = tamfigura;
% hl.Position = [0.6952 0.6683 0.2054 0.1242];

hf.Position = tamfigura3;
hl.Position = [0.7238 0.6671 0.2054 0.0917];



