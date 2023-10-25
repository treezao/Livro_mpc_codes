clear all, 
close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');

Ts = 1; % Periodo de amostragem do processo em segundos

G = 10/(10*s+1); % Modelo por função de transferência nominal do processo
Gz = c2d(G,Ts,'zoh'); % modelo discretizado

Gr = 1/s; % modelo por função de transferência real do sistema
Grz = c2d(Gr,Ts,'zoh');
num = Grz.num{1}; % numerador do modelo discreto
den = Grz.den{1}; % denominador do modelo discreto
na = size(den,2)-1; % ordem do denominador
nb = size(num,2)-2; % ordem do numerador
dd = Grz.inputdelay; % atraso discreto

Nss=80;

Gcoef = step(G,Ts:Ts:Nss*Ts);

%% parâmetros de ajuste

N1 = 1; %horizonte de predição inicial
N2 = 30; % horizonte de predição final
Nu = 5; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle


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
nit = 100 + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidasi = 0*ones(nit,1); % vetor de saídas intermediárias
saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+50:end) = 5;

refs = 0*ones(nit,1); % vetor de referências
refs(nin+10:nit) = 10;

erro = zeros(nit,1); % vetor de erros

%% simulação com referencia futura
for i = nin:nit
    %% simulação do modelo processo
    saidasi(i) = -den(2:end)*saidasi(i-1:-1:i-na) + num*(entradas(i-dd:-1:i-nb-dd-1) );
    saidas(i) = saidasi(i) + perts(i);
    erro(i) = refs(i)-saidas(i);
    
    %% Controlador
    
    Ylivre = Ylivre + Gcoef(1:Nss,1)*duAnt;
    
    eta = saidas(i)-Ylivre(1); %%% erro de predição
    
    Ylivre(1:end-1) = Ylivre(2:end); %%% deslocamento do vetor
    
    %%% resposta livre
    f = Ylivre(N1:N2,1)+eta;
    
    %%% referências
    R = ones(N2-N1+1,1)*refs(i);
    
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
t = ((nin:nit)-nin)*Ts;
vx = (nin:nit);
cores = gray(5);
cores = cores(1:end-1,:);

cores = gray(3);
cores = cores(1:end-1,:);


hf = figure
h=subplot(2,1,1);
plot(30+saidas(vx),'-','LineWidth',tamlinha,'Color',cores(1,:));
grid on
hold
%plot(t(1:nm),y1(1:nm),'-');
plot(30+refs(vx),'--','LineWidth',tamlinha,'Color',cores(2,:));


ylabel('Controlada (%)', 'FontSize', tamletra);
axis([0 100 29 46]);
set(h, 'FontSize', tamletra);
%title('(b)', 'FontSize', tamletra);
%axis([0 pontos -0.1 0.7])
legend('Controlada','Referência','Location','SouthEast')


h=subplot(2,1,2);
plot(47+ entradas(vx),'-','LineWidth',tamlinha,'Color',cores(1,:));
hold
%plot(t(1:nm),u1(1:nm),'-');
grid on
ylabel('Manipulada (%)', 'FontSize', tamletra);
xlabel('Tempo (minutos)', 'FontSize', tamletra);
%axis([0 pontos -2 1])
axis([0 100 44 54]);
set(h, 'FontSize', tamletra);


hf.Position = tamfigura;


