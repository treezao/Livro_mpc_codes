% clear all, 
% close all, 
clc

addpath('../../../Bibliotecas/')
run('../../../Bibliotecas/parametrosFiguras.m')
%%
Ts = 1; % período de amostragem
gz = c2d(tf(1,[2 2 1]),.2); % modelo discreto


A = gz.den{1}; % denominador do modelo do processo
B = gz.num{1}(2:end); % numerador do modelo do processo
Bq = gz.num{1}; % numerador do modelo da perturbação do processo

d=0; % atraso em relação à manipulada
dq = 0; % atraso em relação à perturbação

nb = size(B,2)-1; % ordem do polinômio B
na = size(A,2)-1; % ordem do polinômio A
nbq = size(Bq,2)-1; % ordem do polinômio Bq



%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 15; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 10; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

af = 0; % polo do filtro de referência

umax = 1.2; % valor de controle máximo
umin = -1.2; % valor de controle mínimo

%% montagem das matrizes
Btil = conv(B,[zeros(1,d) 1]); % incorporação do atraso no numerador B

G = zeros(N,N); % matriz dinâmica G0
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
Qe = delta*eye(N);  % matriz de ponderaão dos erros futuros

Kmpc = (G'*Qe*G+Qu)\G'*Qe;

Kmpc1 = Kmpc(1,:);

%% obtenção do controlador equivalente
delta = tf([1 -1],[1 0],-1)
den1 = tf([1 Kmpc1*H],[1 zeros(1,size(Kmpc1*H,2))],-1)
C1num = tf(sum(Kmpc1),1,-1)
C2num = tf(Kmpc1*F,[1 zeros(1,size(Kmpc1*F,2)-1)],-1)

%%% controlador de realimentação
Cr = minreal(C2num/den1/delta)

%%% filtro de referência
Fr = minreal(C1num/C2num)


%%% Função de transferência de malha fechada com e sem filtro de referÊncia
Hr = feedback(Cr*gz,1)

Hrf = minreal(Hr*Fr)

%%% função de transferência de malha fechada para a perturbação
Hq = minreal(feedback(1,Cr*gz)*gz)

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

reffilt = zeros(nit,1);
nfr = size(Fr.den{1},2)-1;

ncr = size(Cr.den{1},2)-1;

%% simulação sem filtro de referência
for k = nin:nit
    %% modelo processo, não mexer
    saidas(k) = -A(2:end)*saidas(k-1:-1:k-na) ...
                  +B*sat(entradas(k-d-1:-1:k-1-nb-d),umin,umax) ...
                  +Bq*perts(k-dq:-1:k-nbq-dq);    
    
    %% -- Controlador PID 
    reffilt(k) = -Fr.den{1}(2:end)*reffilt(k-1:-1:k-nfr) ...
                      +Fr.num{1}*refs(k:-1:k-nfr);

    
    %%% cálculo do sinal de controle do PID
    entradas(k) = -Cr.den{1}(2:end)*entradas(k-1:-1:k-ncr)...
                       +Cr.num{1}*(reffilt(k:-1:k-ncr)-saidas(k:-1:k-ncr));
                  
                  
                  
    %%% cálculo do incremento de controle ótimo
    du(k) = entradas(k)-entradas(k-1);
    
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
hl = legend('PID s/ AW','Referência','Location','NorthEast')
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




