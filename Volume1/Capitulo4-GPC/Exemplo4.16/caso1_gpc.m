% clear all, 
% close all, 
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%
Ts = 1; % período de amostragem
A = poly([0.9512 0.9512]);; % denominador do modelo do processo
B = [0.002418 0.002339]; % numerador do modelo do processo
Bq = [0 B]; % numerador do modelo da perturbação do processo

d=5; % atraso em relação à manipulada
dq = d; % atraso em relação à perturbação

nb = size(B,2)-1; % ordem do polinômio B
na = size(A,2)-1; % ordem do polinômio A
nbq = size(Bq,2)-1; % ordem do polinômio Bq



%% parâmetros de ajuste

N1 = 6; %horizonte de predição inicial
N2 = 50; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 10; % horizonte de controle
Nq = 0; % horizonte da ação antecipativa

delta = 1; % ponderação do erro futuro
lambda = 0.001; % ponderação do esforço de controle

af = 0; % polo do filtro de referência

%% montagem das matrizes
Btil = conv(B,[zeros(1,d) 1]); % incorporação do atraso no numerador B
Bqtil = conv(Bq,[zeros(1,dq) 1]); % incorporação do atraso no numerador Bq


G = zeros(N,N); % matriz dinâmica G0
H = zeros(N,nb+d); % matriz dos incrementos passados de controle 

[E,F] = diofantina(conv(A,[1 -1]),N1,N2); % obtenção dos polinômios Ej, Fj

for i=N1:N2
    EjB = conv(E(i-N1+1,1:i),Btil);
    EjBq = conv(E(i-N1+1,1:i),Bqtil);
    
    G(i-N1+1,1:i) = EjB(i:-1:1);    
    H(i-N1+1,:) = EjB(i+1:end);

end
G = G(:,1:Nu);


G,F,H

%% obtenção do ganho irrestrito
Qu = lambda*eye(Nu); % matriz de ponderação dos incrementos de controle
Qe = delta*eye(N);  % matriz de ponderaão dos erros futuros

Kmpc = (G'*Qe*G+Qu)\G'*Qe;

Kmpc1 = Kmpc(1,:);

%% inicialização vetores
nin = 10;
nit = 100 + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+50:end) = 2;

refs = 0*ones(nit,1); % vetor de referências
refs(nin+10:end) = 1;



rfant = 0;

%% simulação sem filtro de referência
for k = nin:nit
    %% modelo processo, não mexer
    saidas(k) = -A(2:end)*saidas(k-1:-1:k-na) ...
                  +B*entradas(k-d-1:-1:k-1-nb-d) ...
                  +Bq*perts(k-dq:-1:k-nbq-dq);    
    
    %% -- Controlador GPC 
    %%% referencias
    rf = af*rfant + (1-af)*refs(k);
    rfant = rf;
    R = rf*ones(N,1);
    
    %%% cálculo da resposta livre;
    f = F*saidas(k:-1:k-na);
    
    if(~isempty(H))
        f = f+ H*du(k-1:-1:k-nb-d); % parcela dos incrementos de controle
    end
    
    %% Resolve o problema de otimização
    du(k) = Kmpc1*(R-f);
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
ylim([0 1.5])
% h.YTick = [0 0.5 1 1.5];
hl = legend('S/ Filtro','Referência','Location','NorthEast')
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
xlabel('Tempo (amostras)','FontSize', tamletra)
ylabel('\Delta u','FontSize', tamletra)
set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.6952 0.6683 0.2054 0.1242];




