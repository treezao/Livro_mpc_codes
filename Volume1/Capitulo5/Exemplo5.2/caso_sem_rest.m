% clear all, 
% close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');
Ts = 1; % perído de amostragem em minutos

z = tf('z',Ts);
Ts = 0.1;
G = 4/(s+1)^2; % modelo em função de transferência em tempo contínuo
Gz = c2d(G,Ts,'zoh'); % modelo discreto do sistema
num = Gz.num{1}; % numerador do modelo discreto
den = Gz.den{1}; % denominador do modelo discreto
na = size(den,2)-1; % ordem do denominador
nb = size(num,2)-2; % ordem do numerador
dd = Gz.inputdelay; % atraso discreto


%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 90; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 25; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

Nss=100; % horizonte de modelo
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

Ylivre = ones(Nss,1)*0; % 0 é o valor inicial da saída do sistema

%% inicialização vetores
nin = 5;
nit = round(18/Ts) + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+round(12/Ts):end) = -1;

refs = 0*ones(nit,1); % vetor de referências
refs(nin+round(2/Ts):nit) = 5;
refs(nin+round(8/Ts):nit) = 2;


%% simulação sem filtro de referência
for k = nin:nit
    %% modelo processo, não mexer
    saidas(k) = -den(2:end)*saidas(k-1:-1:k-na) + num*(entradas(k-dd:-1:k-nb-dd-1) + perts(k-dd:-1:k-dd-nb-1));
    
    
    %% -- Controlador GDMC 
    Ylivre = Ylivre + Gcoef(1:Nss)*du(k-1);
    
    eta = saidas(k)-Ylivre(1); %%% erro de predição
    
    Ylivre(1:end-1) = Ylivre(2:end); %%% deslocamento do vetor
    
    %%% resposta livre
    f = Ylivre(N1:N2,1)+eta;

    %%% referencias
    R = refs(k)*ones(N,1);
    
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
plot(t,saidas(vx)+52,'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx)+52,'--','LineWidth',tamlinha,'Color',cores(2,:))

hl = legend('S/ restrição','Referência')
ylabel('Nível (%)', 'FontSize', tamletra);
ylim([50 60]);
set(h, 'FontSize', tamletra);
grid on

h = subplot(2,1,2)
plot(t,entradas(vx)+30,'LineWidth',tamlinha,'Color',cores(1,:))
ylabel('Manipulada (%)', 'FontSize', tamletra);
xlabel('Tempo (minutos)', 'FontSize', tamletra);
ylim([20 50]);

grid on

set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7077 0.5411 0.2054 0.1242];




