% clear all, 
% close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');
Ts = 1; % perído de amostragem em minutos

z = tf('z',Ts);
Ts = 1;
Gz = tf([0 0.2 0.4 0.8 0.4 0.2],[1],Ts,'Variable','z^-1'); % modelo discreto do sistema
num = Gz.num{1}; % numerador do modelo discreto
den = Gz.den{1}; % denominador do modelo discreto
na = size(den,2)-1; % ordem do denominador
nb = size(num,2)-2; % ordem do numerador
dd = Gz.inputdelay; % atraso discreto


%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 3; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 2; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = .5; % ponderação do esforço de controle

Nf = 40; % horizonte de modelo filtrado
Nss=5; % horizonte de modelo
Gcoef = step(Gz,Ts:Ts:2*Nss*Ts); % coeficientes da resposta ao degrau
Gqcoef = step(Gz,0:Ts:2*Nss*Ts); % coeficientes da resposta ao degrau da perturbação


%% Monta as matrizes offline que definem o problema DMC
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

%%% calcula a matriz H para o cálculo da resposta livre no caso DMC
H1dmc = [];
H2dmc = [];

H1qdmc = [];
H2qdmc = [];

for i=N1(1):N2(1)
    H1dmc = [H1dmc;Gcoef(i+1:i+Nss)'];
    H2dmc = [H2dmc;Gcoef(1:Nss)'];

    H1qdmc = [H1qdmc;Gqcoef(i+1:i+Nss+1)'];
    H2qdmc = [H2qdmc;Gqcoef(1:Nss+1)'];    
    
end
H = H1dmc-H2dmc
Hq = H1qdmc-H2qdmc

%% cálculo do controlador equivalente
f0 = sum(Kdmc1);
Hu = conv([1 -1],[1 Kdmc1*H])
Huq = Kdmc1*Hq;

Cr = tf(f0,Hu,Ts,'Variable','z^-1') % controlador
Cff = -tf(Huq,[1 Kdmc1*H],Ts,'Variable','z^-1') % controlador antecipativo

mf = feedback(Cr*Gz,1)


%% inicialização vetores
nin = max(Nss,Nf)+1;
nit = 40 + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+round(25/Ts):end) = 2;

refs = 0*ones(nit,1); % vetor de referências
refs(nin+round(4/Ts):end) = 1;


erro = zeros(nit,1); % vetor de erros
yfilt = zeros(nit,N(1)); % vetor com as saidas filtras
dq = zeros(nit,1);

%% simulação sem filtro de referência
for k = nin:nit
    %% modelo processo, não mexer
    saidas(k) = -den(2:end)*saidas(k-1:-1:k-na) + num*(entradas(k-dd:-1:k-nb-dd-1) + perts(k-dd:-1:k-dd-nb-1));
    
    erro(k) = refs(k)-saidas(k);
    
    %% -- Controlador GDMC 
    %%% referencias
    dq(k) = perts(k)-perts(k-1);
    R = refs(k)*ones(N,1);
    f = H*du(k-1:-1:k-Nss) +Hq*dq(k:-1:k-Nss) + saidas(k);
    
    %% Resolve o problema de otimização
    du(k) = Kdmc1*(R-f);
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
% ylim([0 1.6])
% h.YTick = [0 0.5 1 1.5];
hl = legend('DMC','Referência','Location','NorthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
% h.YTick = [-2 -1 0 1 2]
% ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7202 0.4960 0.2054 0.1242]




