% clear all, 
% close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');
Ts = 1; % perído de amostragem em minutos

z = tf('z',Ts);
Gz = 0.1/(z-1.1)/z^2; % modelo discretizado
num = Gz.num{1}; % numerador do modelo discreto
den = Gz.den{1}; % denominador do modelo discreto
na = size(den,2)-1; % ordem do denominador
nb = size(num,2)-2; % ordem do numerador
dd = Gz.inputdelay; % atraso discreto


%% parâmetros de ajuste

N1 = 3; %horizonte de predição inicial
N2 = 25; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 5; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

Nf = 40; % horizonte de modelo filtrado
Nss=80; % horizonte de modelo
Gcoef = step(Gz,Ts:Ts:2*Nss*Ts);

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

H1 = [];
H2 = [];

for i=N1(1):N2(1)
    H1 = [H1;Gcoef(i+1:i+Nss)'];
    H2 = [H2;Gcoef(1:Nss)'];
    
end
H = H1-H2

%% inicialização vetores
nin = max(Nss,Nf)+1;
nit = 150 + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+round(100/Ts):end) = 0.5;

refs = 0*ones(nit,1); % vetor de referências
refs(nin+round(4/Ts):end) = 1;


erro = zeros(nit,1); % vetor de erros
yfilt = zeros(nit,N(1)); % vetor com as saidas filtras


%% simulação sem filtro de referência
for k = nin:nit
    %% modelo processo, não mexer
    saidas(k) = -den(2:end)*saidas(k-1:-1:k-na) + num*(entradas(k-dd:-1:k-nb-dd-1) + perts(k-dd:-1:k-dd-nb-1));
    
    erro(k) = refs(k)-saidas(k);
    
    %% -- Controlador GDMC 
    %%% referencias
    R = refs(k)*ones(N,1);
    f = H*du(k-1:-1:k-Nss) + saidas(k);
    
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
ylim([0 1.6])
h.YTick = [0 0.5 1 1.5];
hl = legend('DMC','Referência','Location','NorthEast')
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
h.YTick = [-2 -1 0 1 2]
ylim([-2.5 2])

ylabel('Manipulada','FontSize', tamletra)
grid on
xlabel('Tempo (amostras)','FontSize', tamletra)

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7202 0.4960 0.2054 0.1242]




