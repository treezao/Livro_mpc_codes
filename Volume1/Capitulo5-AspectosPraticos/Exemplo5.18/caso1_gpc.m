% clear all, 
% close all, 
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');
G = 4/(s*(s+1)^2); % modelo contínuo

Ts = .1; % período de amostragem
Gz = c2d(G,Ts,'zoh');


A = Gz.den{1}; % denominador do modelo do processo
B = Gz.num{1}(2:end); % numerador do modelo do processo
Bq = Gz.num{1}; % numerador do modelo da perturbação do processo

d=0; % atraso em relação à manipulada
dq = 0; % atraso em relação à perturbação

nb = size(B,2)-1; % ordem do polinômio B
na = size(A,2)-1; % ordem do polinômio A
nbq = size(Bq,2)-1; % ordem do polinômio Bq

INTERVALHO_RT_FALHA = 3; %Intervalo de falta de tempo para a otimização 

%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 90; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 25; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = .01; % ponderação do esforço de controle

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

Rbar = G;

%% inicialização vetores
nin = 10;
nit = round(7/Ts) + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
% perts(nin+50:end) = 0.5;

refs = 0*ones(nit+N2,1); % vetor de referências
refs(nin+round(1/Ts):end) = 1;

cont = 0;

%% simulação sem filtro de referência
for k = nin:nit
    %% modelo processo, não mexer
    saidas(k) = -A(2:end)*saidas(k-1:-1:k-na) ...
                  +B*entradas(k-d-1:-1:k-1-nb-d) ...
                  +Bq*perts(k-dq:-1:k-nbq-dq);    
    
    %% -- Controlador GPC 
    %%% referencias
    R = refs(k)*ones(N,1);
    
    %%% cálculo da resposta livre;
    f = F*saidas(k:-1:k-na);
    
    if(~isempty(H))
        f = f+ H*du(k-1:-1:k-nb-d); % parcela dos incrementos de controle
    end
    
    %% Resolve o problema de otimização
    cont = cont + 1;
    if cont == INTERVALHO_RT_FALHA 
        du(k) = 0;
        cont = 0;
    else
        fqp = fqp1*(R-f);

        if(refs(k)==saidas(k))
            [X,FVAL,FLAG(k)] = quadprog(Hqp,fqp);
        elseif(refs(k)>=saidas(k))
            [X,FVAL,FLAG(k)] = quadprog(Hqp,fqp, Rbar,R-f);
        else
            [X,FVAL,FLAG(k)] = quadprog(Hqp,fqp,-Rbar,f-R);
        end
        du(k) = X(1);
    end
    
    entradas(k) = entradas(k-1)+du(k);
    
end

%% plots
t = ((nin:nit)-nin)*Ts;
vx = nin:nit;

cores = gray(3);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(t,saidas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx),'--','LineWidth',tamlinha,'Color',cores(2,:))
% ylim([0 1.5])
% h.YTick = [0 0.5 1 1.5];
hl = legend('GPC 1','Referência','Location','NorthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
% h.YTick = [-2 -1 0 1 2]
% ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on

xlabel('Tempo (segundos)','FontSize', tamletra)

set(h, 'FontSize', tamletra);

hf.Position = tamfigura;
hl.Position = [0.6952 0.6683 0.2054 0.1242];

% hf.Position = tamfigura3;
% hl.Position = [0.7238 0.6671 0.2054 0.0917];



