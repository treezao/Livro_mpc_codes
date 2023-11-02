% clear all, 
% close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');

Ts = 0.05; % Periodo de amostragem do processo em minutos

G = 4/(s+2)^2*exp(-10*Ts*s); % Modelo por função de transferência do processo



Gz = c2d(G,Ts,'zoh'); % modelo discretizado
num = Gz.num{1}; % numerador do modelo discreto
den = Gz.den{1}; % denominador do modelo discreto
na = size(den,2)-1; % ordem do denominador
nb = size(num,2)-2; % ordem do numerador
dd = Gz.inputdelay; % atraso discreto



%% parâmetros de ajuste

N1 = 11; %horizonte de predição inicial
N2 = 40; % horizonte de predição final
Nu = 15; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

umax = 100; % valor máximo do sinal de controle
umin = 0; % valor mínimo do sinal de controle
dumax = 2; % valor máximo do incremento de controle
dumin = -dumax; % valor mínimo do incremento de controle

Nss=80; % horizonte de modelo

Gcoef = step(G,Ts:Ts:Nss*Ts); % coeficientes da resposta ao degrau


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

%% matrizes para o caso com restrições
Hqp = 2*(G'*Qy*G+Qu);
fqp1 = -2*G'*Qy; 

LB = repelem(dumin,Nu')';
UB = repelem(dumax,Nu')';
Rbar = tril(ones(Nu));
Rbar = [Rbar;-Rbar];

%% inicialização vetores
duAnt = 0;% incremento de controle passado
uAnt = 0; % sinal de controle passado

nin = Nss+1;
nit = round(20/Ts) + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit,1); % vetor com as perturbações do sistema
perts(nin+round(1.5/Ts):end) = 20*sin(3*Ts*(nin+round(1.5/Ts):nit));

refs = 0*ones(nit,1); % vetor de referências
refs(nin+round(0.5/Ts):nit) = 50;

erro = zeros(nit,1); % vetor de erros

%% simulação com referencia futura
for i = nin:nit
    %% simulação do modelo processo
    saidas(i) = -den(2:end)*saidas(i-1:-1:i-na) + num*(entradas(i-dd:-1:i-nb-dd-1) + perts(i-dd:-1:i-dd-nb-1));
    
    erro(i) = refs(i)-saidas(i);
    
    %% Controlador
    
    Ylivre = Ylivre + Gcoef(1:Nss,1)*duAnt;
    
    eta = saidas(i)-Ylivre(1); %%% erro de predição
    
    Ylivre(1:end-1) = Ylivre(2:end); %%% deslocamento do vetor
    
    %%% resposta livre
    f = Ylivre(N1:N2,1)+eta;
    
    %%% referências futuras
    R = refs(i);
    
    %%% calculo do incremento de controle ótimo    
    fqp = fqp1*(R-f);
    
    rbar = [repelem(umax-entradas(i-1),Nu)';
             repelem(entradas(i-1)-umin,Nu)';
             ];
    X = quadprog(Hqp,fqp,Rbar,rbar,[],[],LB,UB);
    
    duAtual = X(1);
    
    %%% calculo da ação de controle real
    uAtual = duAtual + uAnt;
    
    %%% aplicar no sistema
    entradas(i) = uAtual;
    du(i) = duAtual;
    
    %%% atualização
    
    uAnt = uAtual;
    duAnt = duAtual;
    
    
end

%% plots
t = ((nin:nit)-nin)*Ts;
vx = nin:nit;


cores = gray(3);
cores = cores(1:end-1,:);

%%% plot dos resultados
hf = figure
h=subplot(3,1,1)
plot(t,saidas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx),'--','Color',cores(1,:),'LineWidth',tamlinha,'Color',cores(2,:))

hl = legend('S/ Bloc.','Referência','Location','NorthEast')
ylabel(['Temperatura',newline,'(^oC)'], 'FontSize', tamletra);
% axis([0 60 0 350]);
grid on
set(h, 'FontSize', tamletra);

h=subplot(3,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
ylabel(['Potência',newline,'(kW)'], 'FontSize', tamletra);
% axis([0 60 0 450]);
grid on
set(h, 'FontSize', tamletra);


h=subplot(3,1,3)
plot(t,du(vx),'LineWidth',tamlinha,'Color',cores(1,:))

ylabel(['\Delta u'], 'FontSize', tamletra);
ylim([-2.5 2.5])
xlabel('Tempo (minutos)', 'FontSize', tamletra);
% axis([0 60 0 450]);
grid on
set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7238 0.6669 0.2054 0.1242]




