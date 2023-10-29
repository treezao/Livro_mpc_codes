clear all, 
close all, 
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%%
s = tf('s');

Ts = 2; % Periodo de amostragem do processo em segundos
z = tf('z',Ts);

G = -5/(50*s+1)*exp(-10*s); % Modelo por função de transferência do processo
Gd = -5/(50*s+1); % modelo sem atraso, este será incorporado depois

Gz = c2d(Gd,Ts,'zoh')*z^-5; % modelo discretizado com atraso


%%% modelos das perturbações
Gq(1,1) = 10/(50*s+1);
Gq(1,2) = 1/(50*s+1);

Gqz = c2d(Gq,Ts,'zoh');


%%% obtenção do modelo em espaço de estados para facilitar simulação
sis = ss([Gz,Gqz]);
nx = size(sis.A,2); % número de estados


%% parâmetros de ajuste
N1 = 6; %horizonte de predição inicial
N2 = 26; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 5; % horizonte de controle
Nq = [5 5]; % horizonte de predição das perturbações

delta = 1/N; % ponderação do erro futuro
lambda = 1/Nu; % ponderação do esforço de controle

umax = [10]; % valor máximo do sinal de controle
umin = [-90]; % valor mínimo do sinal de controle
dumax = [1]; % valor máximo do incremento de controle
dumin = -dumax; % valor mínimo do incremento de controle


Nss=100; % horizonte de modelo

Gcoef = step(G,Ts:Ts:Nss*Ts); % coeficientes da resposta ao degrau
Gqcoef1 = step(Gq(1,1),0:Ts:(Nss-1)*Ts); % coeficientes da resposta ao degrau da perturbação 1
Gqcoef2 = step(Gq(1,2),0:Ts:(Nss-1)*Ts); % coeficientes da resposta ao degrau da perturbação 2


%% montando as matrizes do DMC recursivo

G = zeros(N2,Nu);
G(:,1) = Gcoef(1:N2,1);


for i=2:Nu
    G(i:end,i) = G(1:end-(i-1),1);    
end

G = G(N1:end,:);

Gq1 = zeros(N2,Nq(1));
Gq1(:,1) = Gqcoef1(1:N2,1);

for i=2:Nq(1)
    Gq1(i:end,i) = Gq1(1:end-(i-1),1);    
end

Gq1 = Gq1(N1:end,:);


Gq2 = zeros(N2,Nq(2));
Gq2(:,1) = Gqcoef2(1:N2,1);

for i=2:Nq(2)
    Gq2(i:end,i) = Gq2(1:end-(i-1),1);    
end

Gq2 = Gq2(N1:end,:);


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
nin = 10;
nit = 150 + nin; % número de iterações da simulação

entradas = 0*ones(nit,1); % vetor o sinal de controle
du = zeros(nit,1); % vetor de incrementos de controle

estados = zeros(nx,nit); % vetor de estados do sistema
saidas = 0*ones(nit,1); % vetor com as saídas do sistema

perts = zeros(nit+max(Nq),2); % vetor com as perturbações do sistema
perts(nin+50:end,1) = 5;
perts(nin+100:end,2) = -4;
dq = perts(1:end,:)-[0,0;perts(1:end-1,:)]


refs = 0*ones(nit,1); % vetor de referências
refs(nin+10:end) = 10;

erro = zeros(nit,1); % vetor de erros

%% simulação com referencia futura
for i = nin:nit
    %% simulação do modelo processo
    estados(:,i) = sis.A*estados(:,i-1) + sis.B*[entradas(i-1);perts(i-1,:)'];
    saidas(i) = sis.C*estados(:,i);
    
    erro(i) = refs(i)-saidas(i);
    %% Controlador
    
    Ylivre = Ylivre + Gcoef(1:Nss,1)*du(i-1);
    
    eta = saidas(i)-Ylivre(1); %%% erro de predição
    
    Ylivre(1:end-1) = Ylivre(2:end); %%% deslocamento do vetor
    
    %%% resposta livre
    f = Ylivre(N1:N2,1)+eta;
    
    %%% referências
    R = ones(N2-N1+1,1)*refs(i);
    
    %%% calculo do incremento de controle ótimo    
    % com restrições
    fqp = fqp1*(R-f);
    
    rbar = [repelem(umax-entradas(i-1),Nu)';
             repelem(entradas(i-1)-umin,Nu)'];
    X = quadprog(Hqp,fqp,Rbar,rbar,[],[],LB,UB);
    du(i) = X(1);
    
    % sem restrições
%     du(i) = Kdmc1*(R-f);
    
    %%% aplicar no sistema
    entradas(i) = entradas(i-1)+du(i);    
    
end

%% plots
t = ((nin:nit)-nin)*Ts;
vx = nin:nit;

cores = gray(4);
cores = cores(1:end-1,:);


hf = figure
h=subplot(3,1,1)
plot(t,saidas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
hold on
plot(t,refs(vx),'--','LineWidth',tamlinha,'Color',cores(2,:))
%ylim([0 20])
hl = legend('Controlada','Referência','Location','NorthEast')
% hl.Position = [0.6785 0.7018 0.2368 0.1136];
ylabel('Controlada','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on

h = subplot(3,1,2)
plot(t,entradas(vx),'LineWidth',tamlinha,'Color',cores(1,:))
ylabel('Manipulada','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

h = subplot(3,1,3)
plot(t,du(vx),'LineWidth',tamlinha,'Color',cores(1,:))
ylabel('\Delta u','FontSize', tamletra)
xlabel('Tempo (segundos)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

hf.Position = tamfigura;
hl.Position = [0.6988 0.6618 0.2054 0.1242];

