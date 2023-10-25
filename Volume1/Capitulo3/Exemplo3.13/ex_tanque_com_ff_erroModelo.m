clear all, 
close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');

Ts = 0.1; % Periodo de amostragem do processo em segundos

G = 4/(2*s+1)/(s+1); % Modelo por função de transferência do processo


Gz = c2d(G,Ts,'zoh'); % modelo discretizado
num = Gz.num{1}; % numerador do modelo discreto
den = Gz.den{1}; % denominador do modelo discreto
na = size(den,2)-1; % ordem do denominador
nb = size(num,2)-2; % ordem do numerador
dd = Gz.inputdelay; % atraso discreto


Nss=100;

Gcoef = step(G,Ts:Ts:Nss*Ts);
Gqcoef = step(-G,0:Ts:Nss*Ts);

%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 15; % horizonte de predição final
Nu = 3; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 0.1; % ponderação do esforço de controle

af=0.8; % constante de tempo do filtro de referência

%% montando as matrizes do DMC recursivo

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
duAnt = 0;% incremento de controle passado
uAnt = 0; % sinal de controle passado

nin = Nss+1;
nit = 550 + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+100:end) = -1;
perts(nin+250:end) = 1;

refs = 0*ones(nit,1); % vetor de referências
refs(nin+0:nit) = 5;
refs(nin+250:nit) = 5;
refs(nin+850:end) = 5;

refsf=filter(1-af,[1 -af],refs); % referências filtradas

erro = zeros(nit,1); % vetor de erros

%% simulação com referencia futura
for i = nin:nit-N2-1
    %% simulação do modelo processo
    saidas(i) = -den(2:end)*saidas(i-1:-1:i-na) + num*(entradas(i-dd:-1:i-nb-dd-1) - 0.8* perts(i-dd:-1:i-dd-nb-1));
    
    erro(i) = refs(i)-saidas(i);
    
    %% Controlador
    dq = perts(i)-perts(i-1);
    
    Ylivre = Ylivre + Gcoef(1:Nss,1)*duAnt + Gqcoef(1:Nss,1)*dq;
    
    eta = saidas(i)-Ylivre(1); %%% erro de predição
    
    Ylivre(1:end-1) = Ylivre(2:end); %%% deslocamento do vetor
    
    %%% resposta livre
    f = Ylivre(N1:N2,1)+eta;
    
    %%% referências
    R = ones(N2-N1+1,1)*refsf(i);
    
    %%% calculo do incremento de controle ótimo    
    duAtual = Kdmc1*(R-f);
    
    %%% calculo da ação de controle real
    uAtual = duAtual + uAnt;
    
    %%% aplicar no sistema
    entradas(i) = uAtual;
    
    %%% atualização
    
    uAnt = uAtual;
    duAnt = duAtual;
    
    
end

%% plots
t = (1:nit-N2-1)*Ts;
cores = gray(5);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1)
plot(t,refs(1:nit-N2-1),'--','Color',cores(1,:),'LineWidth',tamlinha)
hold on
plot(t,saidas(1:nit-N2-1),'-','LineWidth',tamlinha,'Color',cores(2,:))
ylabel('Nível (%)', 'FontSize', tamletra);
axis([0 50 0 7]);
grid on

legend('Saída','Referência','Location','SouthEast');

set(h, 'FontSize', tamletra);

h=subplot(2,1,2)
plot(t,entradas(1:nit-N2-1),'-','Color',cores(2,:),'LineWidth',tamlinha)
hold on
plot(t,perts(1:nit-N2-1),'--','LineWidth',tamlinha,'Color',cores(3,:))
ylabel('Entradas (%)', 'FontSize', tamletra);
legend('Controle','Perturbação');
xlabel('Tempo (minutos)', 'FontSize', tamletra);
axis([0 50 -5 10]);
set(h, 'FontSize', tamletra);
grid on

hf.Position = tamfigura;




