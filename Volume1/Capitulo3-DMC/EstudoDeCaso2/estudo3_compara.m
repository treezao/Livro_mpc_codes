%%% Exemplo comparativo entre GDMC e DMC
clear all,
close all,
clc

run('../../../Bibliotecas/parametrosFiguras.m')
%% modelo não-linear e seus parâmetros
global k1 k2 Vc Qf
k1 = 10;
k2 = 10;
Vc = 1;
Qf = 1/30;


Cfbar = 3.288
Cobar = 1.316
%% modelo linearizado
s = tf('s');

Ts = 5; % Periodo de amostragem do processo em segundos
z = tf('z',Ts);

G = 3.433/(103.1*s-1); % Modelo por função de transferência do processo
Gz = c2d(G,Ts,'zoh'); % modelo discretizado

%% parâmetros de ajuste
N1 = 1; %horizonte de predição inicial
N2 = 30; % horizonte de predição final
N = N2-N1+1; % horizonte de predição
Nu = 5; % horizonte de controle

delta = 1; % ponderação do erro futuro
lambda = 1; % ponderação do esforço de controle

umax = [0.05]+Cfbar; % valor máximo do sinal de controle
umin = [-0.05]+Cfbar; % valor mínimo do sinal de controle
dumax = [0.025]; % valor máximo do incremento de controle
dumin = -dumax; % valor mínimo do incremento de controle
ymax = 0.12+Cobar; % valor máximo para a saída
ymin = -0.05+Cobar; % valor mínimo para a saída


Nss=40; % horizonte de modelo
Nf = 25; % horizonte de modelo filtrado

Gcoef = step(G,Ts:Ts:2*Nss*Ts); % coeficientes da resposta ao degrau


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


%% matrizes para o caso com restrições
Hqp = 2*(G'*Qy*G+Qu);
fqp1 = -2*G'*Qy; 

LB = repelem(dumin,Nu')';
UB = repelem(dumax,Nu')';
Rbar = tril(ones(Nu));
Rbar = [Rbar;-Rbar;G;-G];

%% configuração da simulação
nin = Nss+1;
nit = round(400/Ts) + nin; % número de iterações da simulação

perts = zeros(nit,2); % vetor com as perturbações do sistema
perts(nin+round(200/Ts):end) = -.03;

refs = Cobar*ones(nit,1); % vetor de referências
refs(nin+round(50/Ts):end) = 0.1+Cobar;

%% simulações com diferentes betaf
betaf = 0.5; % formato do filtro (az+b)/(z-betaf)^2

estudo3_gdmc

saidas1 = saidas;
entradas1 = entradas;
du1 = du;

%%%
betaf = 0.6; % formato do filtro (az+b)/(z-betaf)^2

estudo3_gdmc

saidas2 = saidas;
entradas2 = entradas;
du2 = du;

%%%
betaf = 0.7; % formato do filtro (az+b)/(z-betaf)^2

estudo3_gdmc

saidas3 = saidas;
entradas3 = entradas;
du3 = du;

%%%
betaf = 0.8; % formato do filtro (az+b)/(z-betaf)^2

estudo3_gdmc

saidas4 = saidas;
entradas4 = entradas;
du4 = du;

%% figuras
t = (0:nit-nin-1)*Ts;

numPlots = 4;
ramp = linspace(0, 0.75, numPlots);
grayColors = [ramp; ramp; ramp]';

hf = figure
h=subplot(2,1,1)
plot(t,saidas1(nin+1:nit),'LineWidth',tamlinha,'Color',grayColors(1,:))
hold on
plot(t,saidas2(nin+1:nit),'LineWidth',tamlinha,'Color',grayColors(2,:))
plot(t,saidas3(nin+1:nit),'LineWidth',tamlinha,'Color',grayColors(3,:))
plot(t,saidas4(nin+1:nit),'LineWidth',tamlinha,'Color',grayColors(4,:))

plot(t,refs(nin+1:nit),'--','LineWidth',tamlinha,'Color','k')
ylim([1.3 1.43])
hl = legend('z_f=0,5','z_f=0,6','z_f=0,7','z_f=0,8','Referência','Location','NorthEast')
ylabel('C_o (mol/L)','FontSize', tamletra)
set(h, 'FontSize', tamletra);
grid on



h = subplot(2,1,2)
plot(t,entradas1(nin+1:nit),'LineWidth',tamlinha,'Color',grayColors(1,:))
hold on
plot(t,entradas2(nin+1:nit),'LineWidth',tamlinha,'Color',grayColors(2,:))
plot(t,entradas3(nin+1:nit),'LineWidth',tamlinha,'Color',grayColors(3,:))
plot(t,entradas4(nin+1:nit),'LineWidth',tamlinha,'Color',grayColors(4,:))


ylabel('C_f (mol/L)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);

xlabel('Tempo (s)','FontSize', tamletra)
grid on
set(h, 'FontSize', tamletra);


hf.Position = tamfigura;
hl.Position = [0.7523 0.4288 0.2054 0.3710];

ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))